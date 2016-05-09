import numpy as np
from test_domain_wall_cobalt import setup_domain_wall_cobalt


def norm(a):
    return np.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def setup_llg_params_near_one(node_count=5, A=3.6 * 4e-7 * np.pi, Ms=6.7e5, K1=4.3, do_precession=True):
    sim = setup_domain_wall_cobalt(
        node_count=node_count, A=A, Ms=Ms, K1=K1, length=1.3, do_precession=do_precession, unit_length=1)
    sim.default_c = 1.23
    sim.gamma = 1.56
    sim.alpha = 2.35
    sim.pins = lambda r: 0
    sim.set_m((0, 0, 1))
    return sim, sim.spin


def compute_jacobean_fd(sim, m, eps=1):
    # use 4th order FD scheme which should produce an exact result without demag
    # hence set eps to 1
    n = sim.mesh.n
    print sim.mesh.coordinates
    # compute the jacobean using the finite difference approximation
    jac = np.zeros((3 * n, 3 * n))
    w = np.array([1. / 12., -2. / 3., 2. / 3., -1. / 12.]) / eps
    once = False
    for j, v in enumerate(np.eye(3 * n)):
        sim.spin[:] = m - 2 * eps * v
        f0 = np.zeros(3 * n)
        sim.sundials_rhs(0, 0, f0)
        sim.spin[:] = m - eps * v
        f1 = np.zeros(3 * n)
        sim.sundials_rhs(0, 0, f1)
        sim.spin[:] = m + eps * v
        f2 = np.zeros(3 * n)
        sim.sundials_rhs(0, 0, f2)
        sim.spin[:] = m + 2 * eps * v
        f3 = np.zeros(3 * n)
        sim.sundials_rhs(0, 0, f3)
        if not once:
            print "m", m
            print "m - 2 * eps * v", m - 2 * eps * v
            print "f0", f0
            print "f1", f1
            print "f2", f2
            print "f3", f3
            once = True
        jac[:, j] = w[0] * f0 + w[1] * f1 + w[2] * f2 + w[3] * f3
    return jac


def compute_jacobean_jtimes(sim, m):
    # use the jtimes function to compute the jacobean
    n = sim.mesh.n
    jac = np.zeros((3 * n, 3 * n))
    tmp = np.zeros(m.shape)
    fy = np.zeros(m.shape)
    jtimes = np.zeros(m.shape)
    for j, v in enumerate(np.eye(3 * n)):
        # use fy=None since it's not used for the computation
        sim.sundials_jtimes(v, jtimes, 0., m, fy)
        jac[:, j] = jtimes
    return jac


def test_compute_fd():
    sim, m = setup_llg_params_near_one()
    # Jacobian computation should be exact with eps=1 or eps=2
    assert np.max(np.abs(compute_jacobean_fd(sim, m, eps=1) - compute_jacobean_fd(sim, m, eps=2))) < 1e-13


def test_compute_jtimes():
    sim, m = setup_llg_params_near_one(node_count=2)
    fd = compute_jacobean_fd(sim, m)
    jtimes = compute_jacobean_jtimes(sim, m)
    print "m"
    print m
    print "FD"
    print fd
    print "JTIMES"
    print jtimes
    assert np.max(np.abs(fd - jtimes)) < 1e-13
