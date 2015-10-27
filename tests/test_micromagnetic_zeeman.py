from fidimag.micro import Zeeman
from fidimag.common import CuboidMesh
from fidimag.micro import Sim
import numpy as np


def varying_field(pos):
    return (1.2 * pos[0], 2.3 * pos[1], 0)


def test_H0_is_indexable_or_callable():
    """
    Test that an exception is raised if H0 is not indexable, and that an
    exception is not raised if H0 is indexable.
    """
    # Test for some different accepted types.
    inputSuccess = ([0., 0., 1.],
                    np.array([0., 0., 1.]),
                    lambda x: x + 0.1)
    for zS in inputSuccess:
        Zeeman(zS)

    # Test for different failing types. Should perhaps use a unittest.TestCase
    # for testing to make this more elegant, but there's probably a reason why
    # it's not used elsewhere.
    inputFailures = [5., -7]
    for zS in inputFailures:
        try:
            Zeeman(zS)
        except ValueError:
            pass
        else:
            raise Exception("Zeeman argument \"{}\" was expected to raise an "
                            "exception, but did not!."
                            .format(zS))

def test_zeeman():

    mesh = CuboidMesh(nx=5, ny=2, nz=1)

    sim = Sim(mesh)
    sim.set_m((1, 0, 0))

    zeeman = Zeeman(varying_field)
    sim.add(zeeman)

    field = zeeman.compute_field()

    assert field[6] == 1.2 * (2 + 0.5)
    assert field[7] == 2.3 * 0.5


if __name__ == "__main__":
    test_zeeman()
    test_H0_is_indexable_or_callable()
