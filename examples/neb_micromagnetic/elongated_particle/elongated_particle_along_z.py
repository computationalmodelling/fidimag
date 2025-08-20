import pytest

"""

Test of the NEBM on an elongated particle of cylindrical shape, based on [1],
where a series of simple magnetic systems are tested using the NEBM in
Spherical coordinates with an Euclidean distance. For the elongated particle
test, the system is a 70 nm long and 12 nm wide (originally it is 13 nm)
nanocylinder, with a uniaxial anisotropy in the direction of the long axis.The
test goal is to find the minimum energy path (MEP) between two degenerate
ferromagnetic states (saturated states along the anisotropy axis), which are
the ground states of this system. The MEP is given by a transverse domain wall
propagation that reverse the ferromagnetic state towards the opposite
anisotropic direction.

In the test, we use the Geodesic and Cartesian NEBM codes (there are issues
with the definition of the angles in the Spherical code). For the Geodesic code
we also describe the magnetisation in Cartesian coordinates but we use a
Geodesic distance as in [2] and a projection of the tangents. In the case of
the Cartesian code we use an Euclidean distance and we do not use projections.
The initial energy band we provide for the algorithm is similar to a coherent
rotation of the spins, which is not the optimal path.

For this particular code, the cylinder long axis (and so the anisotropy axis)
is defined along the z direction. In case we used spherical coordinates,
the angles that define the spin directions are not completely defined when
the spins point in the z direction (poles) thus the algorithm struggles
to find the minimum energy path. It seems this is also the case for Cartesian
coordinates since we haven't projected the tangents.

An interesting test would be to compare the convergences when using a more
symmetrical initial band, like
    initial_images = [[0, 0, -1], [0, 1, 0], [0, 0, 1]]

According to our tests, the Geodesic code is the one that converges faster
and not having much issues with the spin directions.

[1] Dittrich et al., JMMM 250 (2002) L12-L19
[2] Bessarab et al., Computer Physics Communications 196 (2015) 335-347

"""

# FIDIMAG:
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, UniaxialAnisotropy, Demag
from fidimag.common.nebm_geodesic import NEBM_Geodesic
import numpy as np

# Material Parameters ---------------------------------------------------------

A = 10e-12
Kx = 3e5
Ms = 3.98e5

# -----------------------------------------------------------------------------


def relax_neb(sim, k, maxst, simname, initial_images, interpolations,
              save_every=10000, method='Geodesic',
              interpolation_method='linear'
              ):
    """
    Relax a Fidimag simulation using the NEBM:

    sim            :: Simulation object with the system specifications
    k              :: NEBM spring constant
    maxst          :: NEBM max steps for the algorithm evolution
    simname        :: NEBM simulation name
    initial_images :: NEBM images
    interpolations :: A list with the number of inteprolations between
                      the NEBM images
    save_every     :: Save VTK or NPY files every certain number of steps
    method         :: NEBM coordinates: 'Cartesian', 'Geodesic'
    interpolation_method :: Method for interpolating the initial_images,
                            'linear' or 'rotation' (Rodrigues formulae)
    """

    method_dict = {'Geodesic': NEBM_Geodesic}

    neb = method_dict[method](sim,
                              initial_images,
                              interpolations=interpolations,
                              spring_constant=k,
                              name=simname,
                              interpolation_method='rotation',
                              openmp=True
                              )

    neb.relax(max_iterations=maxst,
              save_vtks_every=save_every,
              save_npys_every=save_every,
              stopping_dYdt=1)

    # x, E_interp = neb.compute_polynomial_approximation(200)
    # np.savetxt('interpolated_band_GEODESIC.dat',
    #            np.column_stack((x, E_interp))
    #            )


# Mesh ------------------------------------------------------------------------

# Define an elongated cylinder along the y direction
def cylinder(r, centre, radius):
    if (r[0] - centre[0]) ** 2. + (r[1] - centre[1]) ** 2 <= radius ** 2.:
        return Ms
    else:
        return 0

# Finite differences mesh
mesh = CuboidMesh(nx=6, ny=6, nz=35,
                  dx=2, dy=2, dz=2,
                  unit_length=1e-9)

centre = (np.max(mesh.coordinates[:, 0]) * 0.5,
          np.max(mesh.coordinates[:, 1]) * 0.5)


# Prepare simulation ----------------------------------------------------------

def elongated_part_sim():
    sim = Sim(mesh)
    sim.Ms = lambda r: cylinder(r, centre, 8)
    sim.add(UniformExchange(A=A))
    sim.add(UniaxialAnisotropy(Kx, axis=(0, 0, 1)))  # Anisotropy along y
    sim.add(Demag())

    return sim


# -----------------------------------------------------------------------------

barriers = {}


def test_energy_barrier_cylinder():
    init_im = [(0, 0, -1), (0, 0.9, 0.1), (0, 0, 1)]
    interp = [8, 8]

    for method in ['Geodesic']:
        relax_neb(elongated_part_sim(),
                  1e4, 2000,
                  'neb_cylinder_z-axis_{}'.format(method),
                  init_im,
                  interp,
                  save_every=5000,
                  method=method
                  )

        # Get the energies from the last state
        data = np.loadtxt('neb_cylinder_z-axis_{}_energy.ndt'.format(method))[-1][1:]
        ebarrier = np.abs(np.max(data) - np.min(data)) / (1.602e-19)
        barriers[method] = ebarrier

        # assert ebarrier < 0.017
        # assert ebarrier > 0.005

    for key in barriers.keys():
        print(key, barriers)


if __name__ == '__main__':
    test_energy_barrier_cylinder()
