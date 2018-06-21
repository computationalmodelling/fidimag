# The idea behind these tests is to prevent regression back to a bug which was
# caused by incorrect referencing. This bug caused

import fidimag

def setup_fixture_micro(driver='llg'):
    mesh = fidimag.common.CuboidMesh(nx=10, ny=10, nz=10,
                                     dx=1, dy=1, dz=1,
                                     unit_length=1e-9)
    sim = fidimag.micro.Sim(mesh, driver=driver)
    return mesh, sim

def test_Ms_regression_Exch_micro():
    mesh, sim = setup_fixture_micro()
    A1 = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_m(Ms1)
    exch = fidimag.micro.UniformExchange(A1)
    sim.add(exch)
    sim.set_m(Ms2)
    assert exch.Ms is sim._magnetisation, "The Ms in the Exchange Micro class is not a reference to Ms in the Sim class"
    assert exch.Ms_inv is sim._magnetisation_inv, "The Ms_inv in the Exchange Micro class is not a reference to Ms_inv in the Sim Class"

def test_Ms_regression_Demag_micro():
    mesh, sim = setup_fixture_micro()
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_m(Ms1)
    demag = fidimag.micro.Demag()
    sim.add(demag)
    sim.set_m(Ms2)
    assert demag.Ms is sim._magnetisation, "The Ms in the Demag Micro class is not a reference to Ms in the Sim class"
    assert demag.Ms_inv is sim._magnetisation_inv, "The Ms_inv in the Demag Micro class is not a reference to Ms_inv in the Sim Class"

def test_Ms_regression_SimpleDemag_micro():
    mesh, sim = setup_fixture_micro()
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_m(Ms1)
    demag = fidimag.micro.SimpleDemag()
    sim.add(demag)
    sim.set_m(Ms2)
    assert demag.Ms is sim._magnetisation, "The Ms in the SimpleDemag Micro class is not a reference to Ms in the Sim class"
    assert demag.Ms_inv is sim._magnetisation_inv, "The Ms_inv in the SimpleDemag Micro class is not a reference to Ms_inv in the Sim Class"

def test_Ms_regression_UniaxialAnisotropy_micro():
    mesh, sim = setup_fixture_micro()
    Ku = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_Ms(Ms1)
    anis = fidimag.micro.UniaxialAnisotropy(Ku)
    sim.add(anis)
    sim.set_Ms(Ms2)
    assert anis.Ms is sim._magnetisation, "The Ms in the Exchange Micro class is not a reference to Ms in the Sim class"
    assert anis.Ms_inv is sim._magnetisation_inv, "The Ms_inv in the UniaxialAnisotropy Micro class is not a reference to Ms_inv in the Sim Class"

def test_exchange_Ms_regression_UniaxialAnisotropy():
    mesh, sim = setup_fixture_micro()
    Ku = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_Ms(Ms1)
    anis = fidimag.micro.UniaxialAnisotropy(Ku)
    sim.add(anis)
    sim.set_Ms(Ms2)
    assert anis.Ms is sim._magnetisation, "The Ms in the Exchange Micro class is not a reference to Ms in the Sim class"
    assert anis.Ms_inv is sim._magnetisation_inv, "The Ms_inv in the UniaxialAnisotropy Micro class is not a reference to Ms_inv in the Sim Class"

def test_llg_Ms_regression():
    mesh, sim = setup_fixture_micro(driver='llg')
    Ku = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_Ms(Ms1)
    print(sim.driver._Ms)
    sim.set_Ms(Ms2)
    print(sim.driver._Ms)
    assert sim.driver._Ms is sim._magnetisation

def test_llg_stt_Ms_regression():
    mesh, sim = setup_fixture_micro(driver='llg_stt')
    Ku = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_Ms(Ms1)
    print(sim.driver._Ms)
    sim.set_Ms(Ms2)
    print(sim.driver._Ms)
    assert sim.driver._Ms is sim._magnetisation













def setup_fixture_atomistic(driver='llg'):
    mesh = fidimag.common.CuboidMesh(nx=10, ny=10, nz=10,
                                     dx=1, dy=1, dz=1,
                                     unit_length=1e-9)
    sim = fidimag.atomistic.Sim(mesh, driver=driver)
    return mesh, sim

def test_Ms_regression_Exch_atom():
    mesh, sim = setup_fixture_atomistic()
    A1 = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_mu_s(Ms1)
    exch = fidimag.atomistic.UniformExchange(A1)
    sim.add(exch)
    sim.set_mu_s(Ms2)
    assert exch.mu_s is sim._magnetisation, "The Ms in the Exchange Micro class is not a reference to Ms in the Sim class"
    assert exch.mu_s_inv is sim._magnetisation_inv, "The Ms_inv in the Exchange Micro class is not a reference to Ms_inv in the Sim Class"

def test_Ms_regression_Demag_atom():
    mesh, sim = setup_fixture_atomistic()
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_mu_s(Ms1)
    demag = fidimag.atomistic.Demag()
    sim.add(demag)
    sim.set_mu_s(Ms2)
    assert demag.mu_s is sim._magnetisation, "The Ms in the Demag Micro class is not a reference to Ms in the Sim class"
    assert demag.mu_s_inv is sim._magnetisation_inv, "The Ms_inv in the Demag Micro class is not a reference to Ms_inv in the Sim Class"

def test_Ms_regression_DemagFull_atomistic():
    mesh, sim = setup_fixture_atomistic()
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_mu_s(Ms1)
    demag = fidimag.atomistic.DemagFull()
    sim.add(demag)
    sim.set_mu_s(Ms2)
    assert demag.mu_s is sim._magnetisation, "The Ms in the SimpleDemag Micro class is not a reference to Ms in the Sim class"
    assert demag.mu_s_inv is sim._magnetisation_inv, "The Ms_inv in the SimpleDemag Micro class is not a reference to Ms_inv in the Sim Class"

def test_Ms_regression_Anisotropy_atom():
    mesh, sim = setup_fixture_atomistic()
    Ku = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_mu_s(Ms1)
    anis = fidimag.atomistic.Anisotropy(Ku)
    sim.add(anis)
    sim.set_mu_s(Ms2)
    assert anis.mu_s is sim._magnetisation, "The Ms in the Exchange atomistic class is not a reference to Ms in the Sim class"
    assert anis.mu_s_inv is sim._magnetisation_inv, "The mu_s_inv in the Anisotropy atomistic class is not a reference to mu_s_inv in the Sim Class"

def test_exchange_Ms_regression_CubicAnisotropy_atom():
    mesh, sim = setup_fixture_atomistic()
    Ku = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_mu_s(Ms1)
    anis = fidimag.atomistic.CubicAnisotropy(Ku)
    sim.add(anis)
    sim.set_mu_s(Ms2)
    assert anis.mu_s is sim._magnetisation, "The Ms in the Exchange atomistic class is not a reference to Ms in the Sim class"
    assert anis.mu_s_inv is sim._magnetisation_inv, "The mu_s_inv in the CubicAnisotropy atomistic class is not a reference to Ms_inv in the Sim Class"

def test_llg_mu_s_regression_atom():
    mesh, sim = setup_fixture_atomistic(driver='llg')
    Ku = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_mu_s(Ms1)
    print(sim.driver._mu_s)
    sim.set_mu_s(Ms2)
    print(sim.driver._mu_s)
    assert sim.driver._mu_s is sim._magnetisation

def test_llg_stt_mu_s_regression_atom():
    mesh, sim = setup_fixture_atomistic(driver='llg_stt')
    Ku = 2e-11
    Ms1 = 8.5e5
    Ms2 = 9.0e5
    sim.set_mu_s(Ms1)
    print(sim.driver._mu_s)
    sim.set_mu_s(Ms2)
    print(sim.driver._mu_s)
    assert sim.driver._mu_s is sim._magnetisation
