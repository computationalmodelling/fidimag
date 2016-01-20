import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, Demag, DMI
from fidimag.micro import Zeeman, TimeZeeman
from fidimag.common.fileio import DataReader

mu0 = 4 * np.pi * 1e-7


def compute_field():

    mesh = CuboidMesh(nx=1, ny=1, nz=1, dx=2.0, dy=2.0, dz=2.0, unit_length=1e-9, periodicity=(True, True, False))

    sim = Sim(mesh, name='relax')

    sim.set_tols(rtol=1e-10, atol=1e-14)
    sim.alpha = 0.5
    sim.gamma = 2.211e5
    sim.Ms = 8.6e5
    sim.do_procession = False

    sim.set_m((0,0,1))
    # sim.set_m(np.load('m0.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = Demag(pbc_2d=True)
    sim.add(demag)
    field=demag.compute_field()
    print field

    np.save('m0.npy', sim.spin)


if __name__ == '__main__':

    compute_field()
