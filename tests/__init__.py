from numpy import allclose


def very_close(a, b):
    """ close to machine precision """
    return allclose(a, b, rtol=1e-14, atol=1e-14)
