from libc.math cimport sin, M_PI
import numpy as np 
cimport numpy as np



def fast_sin_init(mesh, double[:] field, *params):
    t, axis, Bmax, fc = params
    for i in range(mesh.n):
        field[3*i+0] = Bmax * axis[0] * sin(fc*t)
        field[3*i+1] = Bmax * axis[1] * sin(fc*t)
        field[3*i+2] = Bmax * axis[2] * sin(fc*t)
