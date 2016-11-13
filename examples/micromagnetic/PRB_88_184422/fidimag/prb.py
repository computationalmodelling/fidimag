import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import dolfin as df
import numpy as np
from micro import Sim
from common import CuboidMesh
from micro import UniformExchange, DMI, UniaxialAnisotropy
from micro import Zeeman, TimeZeeman
from fidimag.common.fileio import DataReader


def test_prb88_184422():
    mu0 = 4 * np.pi * 1e-7

    Ms = 8.6e5
    A = 16e-12
    D = 3.6e-3
    K = 510e3

    mesh = CuboidMesh(nx=100, dx=1, unit_length=1e-9)

    sim = Sim(mesh)

    sim.driver.set_tols(rtol=1e-10, atol=1e-14)

    sim.driver.alpha = 0.5
    sim.driver.gamma = 2.211e5
    sim.Ms = Ms
    sim.do_precession = False

    sim.set_m((0, 0, 1))

    sim.add(UniformExchange(A))
    sim.add(DMI(-D, type='interfacial'))
    sim.add(UniaxialAnisotropy(K, axis=(0, 0, 1)))

    sim.relax(dt=1e-13, stopping_dmdt=0.01, max_steps=5000,
              save_m_steps=None, save_vtk_steps=50)

    m = sim.spin

    mx, my, mz = np.split(m, 3)

    x_array = np.linspace(-49.5, 49.5, 100)

    #plt.plot(x_array, mx)
    #plt.plot(x_array, my)
    #plt.plot(x_array, mz)

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

    #plt.plot(x_array, mx_df)

    assert abs(np.max(mx - mx_df)) < 0.05
    #plt.savefig('dmi_1d.pdf', format='pdf', bbox_inches='tight')

test_prb88_184422()
