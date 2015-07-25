#!/usr/bin/python2

# This script creates multiple fidimag simulation objects and shows how much
# memory (maybe) the simulation objects claim in memory.

import fidimag
import gc
import resource
import sys


def create_simulation():
    return fidimag.micro.Sim(fidimag.micro.FDMesh(x0=-200, dx=1, nx=400,
                                                  y0=-200, dy=1, ny=400,
                                                  unit_length=1e-9, pbc="xy"),
                             name="resource_test")

def measure_memory():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


print("Initial memory before simulation creation: {}MB?".format(measure_memory()/1024.**2))
for zI in range(30):
    sim = create_simulation()
    gc.collect()
    print("Memory after creation of simulation {:02d}: {}MB?"\
          .format(zI + 1, measure_memory()/1024.**2))
    print("\tReference count: {}".format(sys.getrefcount(sim)))
