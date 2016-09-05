from __future__ import print_function
import numpy as np
import pytest
from fidimag.micro.micro_driver import INTEGRATOR_CHOICES
from test_domain_wall_cobalt import setup_domain_wall_cobalt, reference_mz

# inputs to test function are below are of form (integrator, use_jac)
# we use use_jac = False for all
inputs = zip(INTEGRATOR_CHOICES, [False] * len(INTEGRATOR_CHOICES))
# and add a test SUNDIALS with Jacobian as well
inputs.append(("sundials", True))


@pytest.mark.slow
@pytest.mark.parametrize("integrator,use_jac", inputs)
def test_domain_wall_cobalt_integrator_choice(integrator, use_jac):
    sim = setup_domain_wall_cobalt(integrator=integrator, use_jac=use_jac)
    print("Will integrate 1 ns.")
    sim.driver.run_until(1e-9)
    xs, m = sim.mesh.coordinates[:, 0], sim.spin

    # Compare with analytical solution.
    m.shape = (-1, 3)
    mz = m[:, 2]
    ref = [reference_mz(x) for x in xs]
    diff = np.max(np.abs(mz - ref))
    print("Max difference between simulation and reference: {}.".format(diff))

    assert diff < 0.01
