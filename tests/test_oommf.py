import numpy as np
import fidimag.extensions.dipolar as dipolar
from fidimag.micro import FDMesh, UniformExchange, Sim, Demag, DMI
from fidimag.micro.oommf import compute_demag_field, compute_exch_field, compute_dmi_field


def compare_fields(v1, v2):

    v1.shape = (-1, 3)
    v2.shape = (-1, 3)

    f = (v1[:, 0]**2 + v1[:, 1]**2 + v1[:, 2]**2)**0.5

    zero_values = f[:] == 0
    f[zero_values] = 1

    diff = abs(v1 - v2)

    # print 'max error',np.max(diff), np.argmax(diff),len(v1)
    # print v1[0,:]-v2[0,:]

    # print v2[0,:]

    max0 = np.max(diff[:, 0] / f)
    max1 = np.max(diff[:, 1] / f)
    max2 = np.max(diff[:, 2] / f)

    v1.shape = (-1,)
    v2.shape = (-1,)
    return max0, max1, max2


def test_oommf_coefficient():

    res = dipolar.compute_Nxx(10, 1, 1, 1, 2, 3)

    assert abs(-0.000856757528962409449 - res) < 5e-15

    # print clib.compute_Nxx_asy(10,1,1,1,2,3)


def test_exch_field_oommf(A=1e-11, Ms=2.6e5):

    mesh = FDMesh(nx=10, ny=3, nz=2, dx=0.5, unit_length=1e-9)

    sim = Sim(mesh)
    sim.Ms = Ms

    exch = UniformExchange(A=A)
    sim.add(exch)

    def init_m(pos):

        x, y, z = pos

        return (np.sin(x) + y + 2.3 * z, np.cos(x) + y + 1.3 * z, 0)

    sim.set_m(init_m)

    field = exch.compute_field()

    init_m0 = """
    return [list [expr {sin($x*1e9)+$y*1e9+$z*2.3e9}] [expr {cos($x*1e9)+$y*1e9+$z*1.3e9}] 0]
    """
    field_oommf = compute_exch_field(mesh, Ms=Ms, init_m0=init_m0, A=A)

    mx0, mx1, mx2 = compare_fields(field_oommf, field)
    assert max([mx0, mx1, mx2]) < 1e-12


def test_with_oommf_spatial_Ms(A=1e-11):

    def spatial_Ms(pos):
        x, y = pos[0], pos[1]

        if x**2 + y**2 < 5**2:
            return 2e4
        else:
            return 0

    init_m0 = """
    return [list [expr {sin($x*1e9)+$y*1e9+$z*2.3e9}] [expr {cos($x*1e9)+$y*1e9+$z*1.3e9}] 0]
    """

    init_Ms = """

    if { $x*$x + $y*$y < 5e-9*5e-9 } {
        return 2e4
    } else {
        return 0
    }

    """

    mesh = FDMesh(nx=12, ny=10, nz=2, dx=0.5, unit_length=1e-9)

    sim = Sim(mesh)
    sim.Ms = spatial_Ms

    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = Demag()
    sim.add(demag)

    def init_m(pos):

        x, y, z = pos

        return (np.sin(x) + y + 2.3 * z, np.cos(x) + y + 1.3 * z, 0)

    sim.set_m(init_m)

    field = exch.compute_field()
    field_oommf = compute_exch_field(
        mesh, init_m0=init_m0, A=A, spatial_Ms=init_Ms)
    mx0, mx1, mx2 = compare_fields(field_oommf, field)
    assert max([mx0, mx1, mx2]) < 1e-12

    field = demag.compute_field()
    field_oommf = compute_demag_field(
        mesh, spatial_Ms=init_Ms, init_m0=init_m0)

    mx0, mx1, mx2 = compare_fields(field_oommf, field)

    assert max([mx0, mx1, mx2]) < 1e-11


def test_dmi_field_oommf(D=4.1e-3, Ms=2.6e5):

    mesh = FDMesh(nx=10, ny=3, nz=2, dx=0.5, unit_length=1e-9)

    sim = Sim(mesh)
    sim.Ms = Ms

    dmi = DMI(D=D, type='interfacial')
    sim.add(dmi)

    def init_m(pos):

        x, y, z = pos

        return (np.sin(x) + y + 2.3 * z, np.cos(x) + y + 1.3 * z, 0)

    sim.set_m(init_m)

    field = dmi.compute_field()

    init_m0 = """
        return [list [expr {sin($x*1e9)+$y*1e9+$z*2.3e9}] [expr {cos($x*1e9)+$y*1e9+$z*1.3e9}] 0]
        """
    # TODO: check the sign of DMI in OOMMF.
    field_oommf = compute_dmi_field(mesh, Ms=Ms, init_m0=init_m0, D=-D)

    mx0, mx1, mx2 = compare_fields(field_oommf, field)

    assert max([mx0, mx1, mx2]) < 1e-12


def test_demag_field_oommf(Ms=6e5):
    mesh = FDMesh(nx=5, ny=2, nz=3, unit_length=1e-9)
    sim = Sim(mesh)

    sim.Ms = Ms

    demag = Demag()
    sim.add(demag)

    def init_m(pos):

        x = pos[0]

        if x <= 2:
            return (1, 0, 0)
        elif x >= 4:
            return (0, 0, 1)
        else:
            return (0, 1, 0)

    sim.set_m(init_m)
    field = demag.compute_field()
    exact = demag.compute_exact()

    init_m0 = """
    
    if { $x <=2e-9 } {
        return "1 0 0"
    } elseif { $x >= 4e-9 } {
        return "0 0 1"
    } else {
        return "0 1 0"
    }
    """

    field_oommf = compute_demag_field(mesh, Ms=Ms, init_m0=init_m0)

    mx0, mx1, mx2 = compare_fields(field_oommf, exact)
    print mx0, mx1, mx2
    assert max([mx0, mx1, mx2]) < 2e-14

    mx0, mx1, mx2 = compare_fields(field_oommf, field)
    print mx0, mx1, mx2

    assert np.max(abs(field - field_oommf)) < 2e-9


def test_demag_field_oommf_large(Ms=8e5, A=1.3e-11):
    mesh = FDMesh(nx=150, ny=50, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)
    sim = Sim(mesh)

    sim.Ms = Ms

    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = Demag()
    sim.add(demag)

    def init_m(pos):

        x, y, z = pos

        return (np.sin(x) + y + 2.3 * z, np.cos(x) + y + 1.3 * z, 0)

    sim.set_m(init_m)
    demag_field = demag.compute_field()
    exch_field = exch.compute_field()

    #exact = demag.compute_exact()

    init_m0 = """
    return [list [expr {sin($x*1e9)+$y*1e9+$z*2.3e9}] [expr {cos($x*1e9)+$y*1e9+$z*1.3e9}] 0]
    """

    demag_oommf = compute_demag_field(mesh, Ms=Ms, init_m0=init_m0)
    exch_oommf = compute_exch_field(mesh, Ms=Ms, init_m0=init_m0, A=A)

    mx0, mx1, mx2 = compare_fields(demag_oommf, demag_field)
    print mx0, mx1, mx2
    #assert max([mx0,mx1,mx2])< 5e-10

    mx0, mx1, mx2 = compare_fields(exch_oommf, exch_field)
    print mx0, mx1, mx2
    assert max([mx0, mx1, mx2]) < 1e-11

    #mx0,mx1,mx2 = compare_fields(demag_oommf, exact)
    # print mx0,mx1,mx2


def test_energy(Ms=8e5, A=1.3e-11, D=1.32e-3):

    mesh = FDMesh(nx=40, ny=50, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)
    sim = Sim(mesh)

    sim.Ms = Ms

    exch = UniformExchange(A=A)
    sim.add(exch)

    demag = Demag()
    sim.add(demag)

    def init_m(pos):

        x, y, z = pos

        return (np.sin(x) + y + 2.3 * z, np.cos(x) + y + 1.3 * z, 0)

    sim.set_m(init_m)

    demag_energy = demag.compute_energy()
    exch_energy = exch.compute_energy()

    # init_m0="""
    # return [list [expr {sin($x*1e9)+$y*1e9+$z*2.3e9}] [expr {cos($x*1e9)+$y*1e9+$z*1.3e9}] 0]
    #"""

    #field_oommf = compute_exch_field(mesh, Ms=Ms, init_m0=init_m0, A=A)

    exch_energy_oommf = 1.9885853028738599e-19
    demag_energy_oommf = 5.5389695779175673e-19
    dmi_energy_oommf = 2.6657360769014251e-20

    print demag_energy, exch_energy

    assert abs(exch_energy - exch_energy_oommf) / exch_energy_oommf < 3e-15
    assert abs(demag_energy - demag_energy_oommf) / demag_energy_oommf < 1e-10


def test_energy_dmi(Ms=8e5, D=1.32e-3):

    mesh = FDMesh(nx=40, ny=50, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)
    sim = Sim(mesh)

    sim.Ms = Ms

    dmi = DMI(D=D, type='interfacial')
    #dmi = DMI(D=D, type='bulk')
    sim.add(dmi)

    def init_m(pos):

        x, y, z = pos

        return (np.sin(x) + y + 2.3 * z, np.cos(x) + y + 1.3 * z, 1)

    sim.set_m(init_m)

    dmi_energy = dmi.compute_energy()

    # init_m0="""
    # return [list [expr {sin($x*1e9)+$y*1e9+$z*2.3e9}] [expr {cos($x*1e9)+$y*1e9+$z*1.3e9}] 1]
    #"""

    #field_oommf = compute_dmi_field(mesh, Ms=Ms, init_m0=init_m0, D=D)

    dmi_energy_oommf = -4.5665527749090378e-20

    print 'dmi energy', dmi_energy

    assert abs(dmi_energy - dmi_energy_oommf) / dmi_energy_oommf < 1e-15


if __name__ == '__main__':

    # test_demag_field_oommf()

    test_exch_field_oommf()

    # test_demag_field_oommf_large()

    # test_energy()

    # test_energy_dmi()

    # test_dmi_field_oommf()

    #test_with_oommf_spatial_Ms()
