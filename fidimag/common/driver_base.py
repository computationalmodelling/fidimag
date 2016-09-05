from __future__ import division
from __future__ import print_function

import os
import numpy as np


class DriverBase(object):
    """
    Common methods for the micromagnetic and atomistic driver classes
    """

    def __init__(self):
        pass

    # -------------------------------------------------------------------------

    def stat(self):
        return self.integrator.stat()

    def set_default_options(self):
        pass

    def reset_integrator(self, t=0):
        """
        Reset the CVODE integrator and set the simulation time to `t`
        The simulation step is reset to zero
        """
        self.integrator.reset(self.spin, t)
        self.t = t  # also reinitialise the simulation time and step
        self.step = 0

    def set_tols(self, rtol=1e-8, atol=1e-10):
        """
        Set the relative and absolute tolerances for the CVODE integrator
        """
        self.integrator.set_options(rtol, atol)

    def compute_effective_field(self, t):
        """
        Compute the effective field from the simulation interactions,
        calling the method from the corresponding Energy class
        """

        # self.spin[:] = y[:]

        self.field[:] = 0

        for obj in self.interactions:
            self.field += obj.compute_field(t)

    def compute_effective_field_jac(self, t, spin):
        self.field[:] = 0
        for obj in self.interactions:
            if obj.jac:
                self.field += obj.compute_field(t, spin=spin)

    def compute_dmdt(self, dt):
        m0 = self.spin_last
        m1 = self.spin
        dm = (m1 - m0).reshape((3, -1))
        max_dm = np.max(np.sqrt(np.sum(dm ** 2, axis=0)))
        max_dmdt = max_dm / dt
        return max_dmdt

    def run_until(self, t):
        """
        Evolve the system with a micromagnetic driver (LLG, LLG_STT, etc.)
        until a specific time `t`, using the specified integrator.
        The integrator was specified with the right hand side of the
        driver equation

        """

        if t <= self.t:
            if t == self.t and self.t == 0.0:
                self.compute_effective_field(t)
                self.data_saver.save()
            return

        ode = self.integrator

        self.spin_last[:] = self.spin[:]

        flag = ode.run_until(t)

        if flag < 0:
            raise Exception("Run cython run_until failed!!!")

        self.spin[:] = ode.y[:]

        self.t = t
        self.step += 1

        # Update field before saving data
        self.compute_effective_field(t)
        self.data_saver.save()

    def relax(self, dt=1e-11, stopping_dmdt=0.01, max_steps=1000,
              save_m_steps=100, save_vtk_steps=100
              ):
        """

        Evolve the system until meeting the dmdt < stopping_dmdt criteria. We
        can specify to save VTK and NPY files with the magnetisation vector
        field, every certain number of the integrator steps

        TODO: Check what dt is exactly doing

        """

        for i in range(0, max_steps + 1):

            cvode_dt = self.integrator.get_current_step()

            increment_dt = dt

            if cvode_dt > dt:
                increment_dt = cvode_dt

            self.run_until(self.t + increment_dt)

            if save_vtk_steps is not None:
                if i % save_vtk_steps == 0:
                    self.save_vtk()
            if save_m_steps is not None:
                if i % save_m_steps == 0:
                    self.save_m()

            dmdt = self.compute_dmdt(increment_dt)

            print(('step=%d ' +
                   'time=%0.3g ' +
                   'max_dmdt=%0.3g ' +
                   'ode_step=%0.3g') % (self.step,
                                        self.t,
                                        dmdt / self._dmdt_factor,
                                        cvode_dt)
                  )

            if dmdt < stopping_dmdt * self._dmdt_factor:
                break

        if save_m_steps is not None:
            self.save_m()

        if save_vtk_steps is not None:
            self.save_vtk()

    # -------------------------------------------------------------------------
    # Save functions ----------------------------------------------------------
    # -------------------------------------------------------------------------

    def save_vtk(self):
        pass

    def save_m(self):
        """
        Save the magnetisation/spin vector field as a numpy array in
        a NPY file. The files are saved in the `{name}_npys` folder, where
        `{name}` is the simulation name, with the file name `m_{step}.npy`
        where `{step}` is the simulation step (from the integrator)
        """
        if not os.path.exists('%s_npys' % self.name):
            os.makedirs('%s_npys' % self.name)
        name = '%s_npys/m_%g.npy' % (self.name, self.step)
        np.save(name, self.spin)

    def save_skx(self):
        """
        Save the skyrmion number density (sk number per mesh site)
        as a numpy array in a NPY file.
        The files are saved in the `{name}_skx_npys` folder, where
        `{name}` is the simulation name, with the file name `skx_{step}.npy`
        where `{step}` is the simulation step (from the integrator)
        """
        if not os.path.exists('%s_skx_npys' % self.name):
            os.makedirs('%s_skx_npys' % self.name)
        name = '%s_skx_npys/m_%g.npy' % (self.name, self.step)

        # The _skx_number array is defined in the SimBase class in Common
        np.save(name, self._skx_number)
