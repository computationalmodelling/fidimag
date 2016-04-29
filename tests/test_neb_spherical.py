from fidimag.common import neb_spherical
import numpy as np
from numpy.testing import assert_allclose

def test_cartesian2spherical():
    res = neb_spherical.cartesian2spherical(np.array([1,1,0]))
    assert_allclose(res, np.array([np.pi/2, np.pi/4]))

def test_spherical2cartesian():
    # r == 1, so the non-zero coordinates are each sqrt(2)/2
    res = neb_spherical.spherical2cartesian(np.array([np.pi/2, np.pi/4]))
    s = (2**.5)/2
    assert_allclose(res, np.array([s, s, 0]), atol=1e-7)
