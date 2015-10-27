#same as test oommf, the difference is that a working version of oommf is not necessary.
import os
import numpy as np
from fidimag.common import CuboidMesh, UniformExchange, Sim, Demag, DMI
from fidimag.micro.oommf import compute_demag_field, compute_exch_field, compute_dmi_field
from fidimag.micro.omf import OMF2


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



def test_exch_field_oommf(A=1e-11, Ms=2.6e5):

    mesh = CuboidMesh(nx=10, ny=3, nz=2, dx=0.5, unit_length=1e-9)

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
    omf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'omfs',
                            'test_exch_field_oommf.ohf'
                            )
    ovf = OMF2(omf_file)
    field_oommf = ovf.get_all_mags()

    #field_oommf = compute_exch_field(mesh, Ms=Ms, init_m0=init_m0, A=A)

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

    mesh = CuboidMesh(nx=12, ny=10, nz=2, dx=0.5, unit_length=1e-9)

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
    #field_oommf = compute_exch_field(mesh, init_m0=init_m0, A=A, spatial_Ms=init_Ms)

    omf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'omfs','test_with_oommf_spatial_Ms_Exchange.ohf')
    ovf = OMF2(omf_file)
    field_oommf = ovf.get_all_mags()

    mx0, mx1, mx2 = compare_fields(field_oommf, field)
    assert max([mx0, mx1, mx2]) < 1e-12

    field = demag.compute_field()
    #field_oommf = compute_demag_field(mesh, spatial_Ms=init_Ms, init_m0=init_m0)
    omf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'omfs','test_with_oommf_spatial_Ms_Demag.ohf')
    ovf = OMF2(omf_file)
    field_oommf = ovf.get_all_mags()

    mx0, mx1, mx2 = compare_fields(field_oommf, field)

    assert max([mx0, mx1, mx2]) < 1e-11


def test_dmi_field_oommf(D=4.1e-3, Ms=2.6e5):

    mesh = CuboidMesh(nx=10, ny=3, nz=2, dx=0.5, unit_length=1e-9)

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
    #field_oommf = compute_dmi_field(mesh, Ms=Ms, init_m0=init_m0, D=-D)
    omf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'omfs','test_dmi_field_oommf.ohf')
    ovf = OMF2(omf_file)
    field_oommf = ovf.get_all_mags()

    mx0, mx1, mx2 = compare_fields(field_oommf, field)

    assert max([mx0, mx1, mx2]) < 1e-12



def test_demag_field_oommf_large(Ms=8e5, A=1.3e-11):
    mesh = CuboidMesh(nx=150, ny=50, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)
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

    #demag_oommf = compute_demag_field(mesh, Ms=Ms, init_m0=init_m0)
    #exch_oommf = compute_exch_field(mesh, Ms=Ms, init_m0=init_m0, A=A)

    omf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'omfs','test_demag_field_oommf_large_Demag.ohf')
    ovf = OMF2(omf_file)
    demag_oommf = ovf.get_all_mags()

    omf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'omfs','test_demag_field_oommf_large_Exchange.ohf')
    ovf = OMF2(omf_file)
    exch_oommf = ovf.get_all_mags()

    mx0, mx1, mx2 = compare_fields(demag_oommf, demag_field)
    #print mx0, mx1, mx2
    assert max([mx0,mx1,mx2])< 5e-10

    mx0, mx1, mx2 = compare_fields(exch_oommf, exch_field)
    #print mx0, mx1, mx2
    assert max([mx0, mx1, mx2]) < 1e-11



if __name__ == '__main__':

    #test_exch_field_oommf()

    #test_with_oommf_spatial_Ms()

    test_demag_field_oommf_large()


    #test_dmi_field_oommf()

    #test_with_oommf_spatial_Ms()

