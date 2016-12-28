from __future__ import print_function
from fidimag.common import CuboidMesh
from fidimag.micro import Sim
from fidimag.micro import Zeeman
from fidimag.micro import UniaxialAnisotropy
import numpy as np
import unittest

def run_sim():
	mesh = CuboidMesh()
	sim = Sim(mesh, name='spin')
	alpha = 0.1
	gamma = 2.21e5
	sim.alpha = alpha
	sim.driver.gamma = gamma
	sim.mu_s = 1.0

	sim.set_m((1, 0, 0))
	H0 = 1e5
	sim.add(Zeeman((0, 0, H0)))
	sim.driver.run_until(1e-10)
	sim.driver.run_until(0.5e-10)

class mytest(unittest.TestCase):
	def test_raises_valueerror_time_integration(self):
		self.assertRaises(ValueError, run_sim)
	 