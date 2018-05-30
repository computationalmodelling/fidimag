cimport numpy as np
import numpy as np

# -----------------------------------------------------------------------------

cdef extern from "common_clib.h":

    # From: llg.c
    void llg_rhs(double * dm_dt, double * spin,
                 double *h, double *alpha, int *pins,
                 double gamma, int n, int do_precession, double default_c)


    void llg_rhs_jtimes(double *jtn, double *m, double *h,
                        double *mp, double *hp, double *alpha, int *pins,
                        double gamma, int n,
                        int do_precession, double default_c)

    void compute_stt_field_c(double *spin, double *field,
                             double *jx, double *jy, double *jz,
                             double dx, double dy, double dz, int *ngbs, int n)

    # From: stt.c
    void llg_stt_rhs(double *dm_dt, double *m, double *h, double *h_stt,
                     double *alpha,double beta, double u0, double gamma, int n)

    void llg_stt_cpp(double *dm_dt, double *m, double *h, double *p,
                     double *alpha, int *pins, double *a_J,
                     double beta, double gamma, int n)

    # -------------------------------------------------------------------------
    # From steepest_descent.c

    void sd_update_spin (double *spin, double *spin_last, double *magnetisation,
                         double *mxH, double *mxmxH, double *mxmxH_last, double *tau,
                         int* pins, int n)

    void sd_compute_step (double *spin, double *spin_last, double *magnetisation,
                          double *field,
                          double *mxH, double *mxmxH, double *mxmxH_last, double *tau,
                          int *pins, int n, int counter, double tmin, double tmax)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def compute_llg_rhs(double [:] dm_dt,
                    double [:] spin,
                    double [:] field,
                    double [:] alpha,
                    int [:] pins,
                    gamma, n, do_precession, default_c):
    llg_rhs(&dm_dt[0], &spin[0], &field[0], &alpha[0], &pins[0],
            gamma, n, do_precession, default_c)


def compute_llg_jtimes(double [:] jtn,
                       double [:] m,
                       double [:] field,
                       double [:] mp,
                       double [:] field_p,
                       double [:] alpha,
                       int [:] pins,
                       gamma, n, do_precession, default_c):
    llg_rhs_jtimes(&jtn[0], &m[0], &field[0], &mp[0], &field_p[0],
                   &alpha[0], &pins[0], gamma, n, do_precession, default_c)

# -----------------------------------------------------------------------------

def compute_stt_field(double [:] spin,
                      double [:] field,
                      double [:] jx,
                      double [:] jy,
                      double [:] jz,
                      dx, dy, dz,
                      int [:, :] ngbs,
                      n
                      ):
    compute_stt_field_c(&spin[0], &field[0], &jx[0], &jy[0],&jz[0],
                        dx, dy, dz, &ngbs[0, 0], n)

def compute_llg_stt_rhs(double [:] dm_dt,
                        double [:] spin,
                        double [:] field,
                        double [:] field_stt,
                        double [:] alpha,
                        beta, u0, gamma, n):
    llg_stt_rhs(&dm_dt[0], &spin[0], &field[0], &field_stt[0],
                &alpha[0], beta, u0, gamma, n)



def compute_llg_stt_cpp(double [:] dm_dt,
                        double [:] spin,
                        double [:] field,
                        double [:] p,
                        double [:] alpha,
                        int [:] pin,
                        double [:] a_J,
                        beta, gamma, n):
    llg_stt_cpp(&dm_dt[0], &spin[0], &field[0], &p[0],
                &alpha[0], &pin[0], &a_J[0], beta, gamma, n)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def compute_sd_spin(double [:] spin,
                    double [:] spin_last,
                    double [:] magnetisation,
                    double [:] mxH,
                    double [:] mxmxH,
                    double [:] mxmxH_last,
                    double [:] tau,
                    int [:] pins,
                    n):

    sd_update_spin(&spin[0], &spin_last[0], &magnetisation[0], &mxH[0],
                   &mxmxH[0], &mxmxH_last[0], &tau[0], &pins[0], n
                   )

def compute_sd_step(double [:] spin,
                    double [:] spin_last,
                    double [:] magnetisation,
                    double [:] field,
                    double [:] mxH,
                    double [:] mxmxH,
                    double [:] mxmxH_last,
                    double [:] tau,
                    int [:] pins,
                    n, counter, tmin, tmax):

    sd_compute_step(&spin[0], &spin_last[0], &magnetisation[0],
                    &field[0], &mxH[0],
                    &mxmxH[0], &mxmxH_last[0], &tau[0], &pins[0],
                    n, counter, tmin, tmax
                    )

def normalise(a):
    """
    normalise the given array a
    """
    a.shape = (-1, 3)
    b = np.sqrt(a[:, 0] ** 2 + a[:, 1] ** 2 + a[:, 2] ** 2)
    ids = (b == 0)
    b[ids] = 1.0
    a[:, 0] /= b
    a[:, 1] /= b
    a[:, 2] /= b
    a.shape = (-1,)

def init_scalar(value, mesh, *args):

    n = mesh.n

    mesh_v = np.zeros(n)

    if isinstance(value, (int, float)):
        mesh_v[:] = value
    elif hasattr(value, '__call__'):
        for i in range(n):
            mesh_v[i] = value(mesh.coordinates[i], *args)

    elif isinstance(value, np.ndarray):
        if value.shape == mesh_v.shape:
            mesh_v[:] = value[:]
        else:
            raise ValueError("Array size must match the mesh size")

    return mesh_v

def init_vector(m0, mesh, norm=False, *args):

    n = mesh.n

    spin = np.zeros((n, 3))

    if isinstance(m0, list) or isinstance(m0, tuple):
        spin[:, :] = m0
        spin = np.reshape(spin, 3 * n, order='C')

    elif hasattr(m0, '__call__'):
        v = m0(mesh.coordinates[0], *args)
        if len(v) != 3:
            raise Exception(
                'The length of the value in init_vector method must be 3.')
        for i in range(n):
            spin[i, :] = m0(mesh.coordinates[i], *args)
        spin = np.reshape(spin, 3 * n, order='C')

    elif isinstance(m0, np.ndarray):
        if m0.shape == (3, ):
            spin[:] = m0  # broadcasting
        else:
            spin.shape = (-1)
            spin[:] = m0  # overwriting the whole thing

    spin.shape = (-1,)

    if norm:
        normalise(spin)

    return spin

def init_vector_func_fast(m0, mesh, double[:] field, norm=False, *args):
    """
    An unsafe method of setting the field. Depends on
    the setter code being memory safe.

    m0 must be a Python function that takes the mesh
    and field as arguments. Within that, the user
    must handle evaluating the function at different
    coordinate points. It needs to be able to handle
    the spatial dependence itself, and write the
    field valuse into the field array. This can
    be written with Cython which will give much
    better performance. For example:

    from libc.math cimport sin

    def fast_sin_init(mesh, double[:] field, *params):
        t, axis, Bmax, fc = params
        for i in range(mesh.n):
            field[3*i+0] = Bmax * axis[0] * sin(fc*t)
            field[3*i+1] = Bmax * axis[1] * sin(fc*t)
            field[3*i+2] = Bmax * axis[2] * sin(fc*t)

    """  
    m0(mesh, field, *args)
    if norm:
        normalise(field)
    return field
