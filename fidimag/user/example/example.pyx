from libc.math cimport cos, sin

def fast_sin_init(mesh, double[:] field, *params):
    t, axis, Bmax, fc = params
    for i in range(mesh.n):
        field[3*i+0] = Bmax * axis[0] * sin(fc*t)
        field[3*i+1] = Bmax * axis[1] * sin(fc*t)
        field[3*i+2] = Bmax * axis[2] * sin(fc*t)

def TimeZeemanFast_test_time_fun(mesh, double[:] field, *params):
    cdef int i
    t, frequency = params
    for i in range(mesh.n):
        field[3*i+0] = 0
        field[3*i+1] = 0
        field[3*i+2] = 10 * cos(frequency * t)
