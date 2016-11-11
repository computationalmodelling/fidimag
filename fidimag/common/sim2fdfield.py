import finitedifferencefield
from fidimag.common import CuboidMesh

def fidimag_to_finitedifferencefield(sim):
    """
    fidimag_to_finitedifferencefield(sim)

    This function takes a Fidimag simulation object, and constructs a
    Finite Difference Field object which has the magnetisation configuration
    from the simulation at the last time step.
    """
    cmin = np.array([sim.mesh.x0, sim.mesh.y0, sim.mesh.z0])*sim.mesh.unit_length
    cmax = tuple(cmin +np.array([mesh.Lx, mesh.Ly, mesh.Lz])*sim.mesh.unit_length)
    cmin = tuple(cmin)
    d = tuple(np.array([sim.mesh.dx, sim.mesh.dy, sim.mesh.dz])*sim.mesh.unit_length)
    field = finitedifferencefield.Field(cmin, cmax, d)
    numpyfield = sim.spin.copy().reshape((mesh.nx, mesh.ny, mesh.nz, 3))
    field.f = numpyfield
    return field

