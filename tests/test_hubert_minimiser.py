import fidimag
import fidimag.common.constant as C
import numpy as np


def test_hubert_minimiser_1D_DW():

    nx, ny, nz = 100, 1, 1
    dx, dy, dz = 1, 1, 1

    mesh = fidimag.common.CuboidMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz,
                                    periodicity=(False, False, False),
                                    unit_length=1e-9)

    sim = fidimag.micro.Sim(mesh, name='1Dmicro', driver='hubert_minimiser')

    Ms = 0.86e6
    A = 13e-12
    Ku = 0.4e6

    # Define the magnetisation
    sim.set_Ms(Ms)

    # Add the magnetic interactions
    sim.add(fidimag.micro.UniformExchange(A))
    sim.add(fidimag.micro.UniaxialAnisotropy(Ku, axis=(0, 0, 1)))
    # sim.add(fidimag.atomistic.DMI(D, dmi_type='interfacial'))
    # sim.add(fidimag.micro.Zeeman((0, 0, B)))

    xs = mesh.coordinates[:, 0]
    centre_x = (xs.max() + xs.min()) * 0.5 + xs.min()

    def m_initial(r):
        x, y, z = r[0], r[1], r[2]
        if x < centre_x:
            return (0, 0.1, -.9)
        else:
            return (0, 0.1, .9)
    sim.set_m(m_initial)

    mesh_vol = mesh.n * mesh.dx * mesh.dy * mesh.dz * 1e-27
    Kd = C.mu_0 * (Ms ** 2) * 0.5 * mesh_vol
    # Scale the energy for the minimiser
    sim.driver.energyScale = Kd

    sim.driver.minimise(stopping_dE=1e-6, maxCreep=6,
                        eta_scale=1e-6, log_steps=10)
    mz = sim.spin.reshape(-1, 3)[:, 2]

    # ANALYSIS

    def mz_dw_analyt(x, xc):
        """Analytical DW function using DW width param"""
        # DW width from theory
        deltaB = np.sqrt(A / Ku)
        return np.tanh((x - xc) / deltaB)

    # Analytical model in nm units
    x = sim.mesh.coordinates[:, 0]
    y = mz_dw_analyt(x * 1e-9, 50 * 1e-9)

    # Mean absolute error at points within the DW width:
    # The maximum "deviation" a spin might have from the curve is 2, since
    # the spin-z is in the [-1, 1] range
    deltaB = np.sqrt(A / Ku) * 1e9
    xc = 50  # sample centre
    ftr = np.logical_and(x <= xc + deltaB, x >= xc - deltaB)

    MAE_dw = np.mean(np.abs(y[ftr] - mz[ftr]))
    print('Domain wall MAE', MAE_dw)
    # This is at least 1% the max deviation
    assert MAE_dw < 0.02


if __name__ == "__main__":
    test_hubert_minimiser_1D_DW()
