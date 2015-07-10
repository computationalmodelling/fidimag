import fidimag.micro
import numpy as np

xMax = 50.  # Upper bound of the domain, centred about zero. Domain only varies
            # in x.
D = 4e-3  # DMI coefficient, [Jm-2].

def m_init(pos):
    """Define magnetisation texture corresponding to +Mz at the edges, -Mz at
    centre, and rotating through My. Encourages helices that change chirality
    at the centre.."""
    Mz = abs(pos[0]) * 2. / xMax - 1
    My = (1 - Mz ** 2.) ** 0.5
    My *= 1 if pos[0] < 0 else -1  # Change direction! Could be the wrong way
                                   # though...
    return [0., My, Mz]

def dmi_variance(pos):
    """Define how the DMI is defined on the simulation domain. Is basically:

    ----------|----------
         D         -D

    """
    return D if pos[0] < 0 else -D

# Create a mesh that looks like this: ------------------
mesh = fidimag.micro.FDMesh(nx=100, dx=1, x0=-50, unit_length=1e-9)
sim = fidimag.micro.Sim(mesh)

# Dynamics parameters.
sim.set_tols(rtol=1e-6, atol=1e-6)
sim.alpha = 0.5
sim.gamma = 2.211e5
sim.Ms = 8.6e5
sim.do_procession = False

# Set magnetisation and add energies.
sim.set_m(m_init)
sim.add(fidimag.micro.UniformExchange(A=1.3e-11))
sim.add(fidimag.micro.DMI(D=dmi_variance))

# Relax and save.
sim.relax(dt=1e-13, stopping_dmdt=1e-2, save_m_steps=None, save_vtk_steps=50)
np.save('m0.npy', sim.spin)
