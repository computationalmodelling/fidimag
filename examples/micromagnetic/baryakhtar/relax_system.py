import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fidimag.common import CuboidMesh
from fidimag.micro import Sim
from fidimag.micro import Demag
from fidimag.micro import Zeeman
from fidimag.micro import UniformExchange


def relax_system(mesh):

    sim = Sim(mesh, chi=1e-3, name='relax', driver='llbar_full')

    sim.set_tols(rtol=1e-7, atol=1e-7)
    sim.Ms = 8.0e5
    sim.alpha = 0.1
    sim.beta = 0
    sim.gamma = 2.211e5

    sim.set_m((1, 0.25, 0.1))
    # sim.set_m(np.load('m0.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    mT = 795.7747154594767
    zeeman = Zeeman([-100 * mT, 4.3 * mT, 0], name='H')
    sim.add(zeeman, save_field=True)

    demag = Demag()
    sim.add(demag)

    ONE_DEGREE_PER_NS = 17453292.52

    sim.relax(dt=1e-12, stopping_dmdt=0.01,
              max_steps=5000, save_m_steps=100, save_vtk_steps=50)

    np.save('m0.npy', sim.spin)

if __name__ == "__main__":

    mesh = CuboidMesh(nx=20, ny=20, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)
    relax_system(mesh)
