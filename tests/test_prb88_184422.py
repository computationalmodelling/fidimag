import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.micro import Sim
from fidimag.micro import FDMesh
from fidimag.micro import UniformExchange, DMI, UniaxialAnisotropy
from fidimag.micro import Zeeman, TimeZeeman
from fidimag.common.fileio import DataReader
from fidimag.micro.omf import OMF2
import os

# def test_prb88_184422():


def run_fidimag(mesh):

    mu0 = 4 * np.pi * 1e-7

    Ms = 8.6e5
    A = 16e-12
    D = -3.6e-3
    K = 510e3

    sim = Sim(mesh)

    sim.set_tols(rtol=1e-10, atol=1e-10)

    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = Ms
    sim.do_precession = False

    sim.set_m((0, 0, 1))

    sim.add(UniformExchange(A))
    sim.add(DMI(D, type='interfacial'))
    sim.add(UniaxialAnisotropy(K, axis=(0, 0, 1)))

    sim.relax(dt=1e-13, stopping_dmdt=0.01, max_steps=5000,
              save_m_steps=None, save_vtk_steps=50)

    m = sim.spin
    return m.copy()


def run_dolfin():

    import dolfin as df

    x_array = np.linspace(-49.5, 49.5, 100)

    mesh = df.IntervalMesh(100, -50, 50)

    Delta = np.sqrt(A / K)
    xi = 2 * A / D

    Delta_s = Delta * 1e9

    V = df.FunctionSpace(mesh, "Lagrange", 1)
    u = df.TrialFunction(V)
    v = df.TestFunction(V)
    u_ = df.Function(V)
    F = -df.inner(df.nabla_grad(u), df.nabla_grad(v)) * df.dx - \
        (0.5 / Delta_s**2) * df.sin(2 * u) * v * df.dx
    F = df.action(F, u_)

    J = df.derivative(F, u_, u)

    # the boundary condition is from equation (8)
    theta0 = np.arcsin(Delta / xi)
    ss = 'x[0]<0? %g: %g ' % (-theta0, theta0)

    u0 = df.Expression(ss)

    def u0_boundary(x, on_boundary):
        return on_boundary

    bc = df.DirichletBC(V, u0, u0_boundary)

    problem = df.NonlinearVariationalProblem(F, u_, bcs=bc, J=J)
    solver = df.NonlinearVariationalSolver(problem)
    solver.solve()

    u_array = u_.vector().array()

    mx_df = []
    for x in x_array:
        mx_df.append(u_(x))

    return mx_df


def test_prb88_184422():
    mesh = FDMesh(nx=100, dx=1, unit_length=1e-9)
    Ms = 8.6e5
    m = run_fidimag(mesh)
    omf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'omfs',
                            'dmi-Oxs_TimeDriver-Magnetization-00-0000963.omf')

    ovf = OMF2(omf_file)
    m_oommf = ovf.get_all_mags()
    assert max(abs(m_oommf / Ms - m)) < 5e-7

if __name__ == '__main__':
	test_prb88_184422()
