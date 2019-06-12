# distutils: language = c++

cimport numpy as np
import numpy as np

# C++ modules
from libcpp cimport bool

# -----------------------------------------------------------------------------

cdef extern from "c_clib.h":

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
                         double *mxH, double *mxmxH, double *mxmxH_last, double tau,
                         int* pins, int n)

    void sd_compute_step (double *spin, double *spin_last, double *magnetisation,
                          double *field,
                          double *mxH, double *mxmxH, double *mxmxH_last, double tau,
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
                    double tau,
                    int [:] pins,
                    n):

    sd_update_spin(&spin[0], &spin_last[0], &magnetisation[0], &mxH[0],
                   &mxmxH[0], &mxmxH_last[0], tau, &pins[0], n
                   )

def compute_sd_step(double [:] spin,
                    double [:] spin_last,
                    double [:] magnetisation,
                    double [:] field,
                    double [:] mxH,
                    double [:] mxmxH,
                    double [:] mxmxH_last,
                    double tau,
                    int [:] pins,
                    n, counter, tmin, tmax):

    sd_compute_step(&spin[0], &spin_last[0], &magnetisation[0],
                    &field[0], &mxH[0],
                    &mxmxH[0], &mxmxH_last[0], tau, &pins[0],
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

def init_vector(m0, mesh, dim=3, norm=False, *args):
    n = mesh.n
    spin = np.zeros((n, dim))
    if isinstance(m0, list) or isinstance(m0, tuple):
        spin[:, :] = m0
        spin = np.reshape(spin, dim * n, order='C')
    elif hasattr(m0, '__call__'):
        v = m0(mesh.coordinates[0], *args)
        if len(v) != dim:
            raise Exception(
                'The length of the value in init_vector method must be {}.'.format(dim))
        for i in range(n):
            spin[i, :] = m0(mesh.coordinates[i], *args)
        spin = np.reshape(spin, dim * n, order='C')
    elif isinstance(m0, np.ndarray):
        if m0.shape == (dim, ):
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

    An unsafe method of setting the field. Depends on the setter code being
    memory safe.

    m0 must be a Python function that takes the mesh and field as arguments.
    Within that, the user must handle evaluating the function at different
    coordinate points. It needs to be able to handle the spatial dependence
    itself, and write the field valuse into the field array. This can be
    written with Cython which will give much better performance. For example:

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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# CLASSES

# -----------------------------------------------------------------------------
# C++ definitions

# See:
# https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html
#
# - Methods that are redefined in inherited classes are only specified in the
# base class
# - These can be declared in a .pxd file:

cdef extern from "c_energy.h":

    cdef cppclass Energy:
        # except +: Without this declaration, C++ exceptions originating from
        # the constructor will not be handled by Cython.
        Energy() except +

        void compute_field(double t)
        double compute_energy()

        # # Not using these variables from Cython:
        # bool set_up
        # int nx, ny, nz, n
        # double dx, dy, dz
        # double unit_length
        # double *spin
        # double *Ms
        # double *Ms_inv
        # double *field
        # double *energy
        # double *coordinates
        # int *ngbs
        int interaction_id

    cdef cppclass ExchangeEnergy(Energy):
        ExchangeEnergy() except +

        void setup(double * A, MicroSim * sim)
        # double *A

# cdef extern from "c_energy.cpp":
#     pass

# -----------------------------------------------------------------------------
# Python definitions

# Cython initializes C++ class attributes of a cdef class using the nullary
# constructor. If the class you’re wrapping does not have a nullary
# constructor, you must store a pointer to the wrapped class and manually
# allocate and deallocate it. A convenient and safe place to do so is in the
# __cinit__ and __dealloc__ methods which are guaranteed to be called exactly
# once upon creation and deletion of the Python instance.

cdef class PyEnergy:
    cdef Energy *thisptr
    def __cinit__(self):
        # No need to allocate memory, we are not using the C++ base class
        # self.thisptr = new Energy()
        if type(self) != PyEnergy:
            return
        print("In Python A")

    # Necessary?:
    def compute_field(self, t):
        self.thisptr.compute_field(t)

    def compute_energy(self, time):
        return self.thisptr.compute_energy()

    def add_interaction_to_sim(self, PyMicroSim sim):
        sim.thisptr.add_interaction(<void *> self.thisptr,
                                    self.thisptr.interaction_id)

    def get_interaction_id(self):
        return self._thisptr.interaction_id


cdef class PyExchangeEnergy(PyEnergy):
    cdef ExchangeEnergy *_thisptr
    # Try cinit:
    def __cinit__(self, double [:] A, PyMicroSim sim):
        print("In Python B")

        self._thisptr = self.thisptr = new ExchangeEnergy()
        self._thisptr.setup(&A[0], sim.thisptr)

    def setup(self, double [:] A, PyMicroSim sim):
        self._thisptr.setup(&A[0], sim.thisptr)

    def __dealloc__(self):
        del self._thisptr


# Simulation class ------------------------------------------------------------

cdef extern from "c_micro_sim.h":

    cdef cppclass MicroSim:
        # except +: Without this declaration, C++ exceptions originating from
        # the constructor will not be handled by Cython.
        MicroSim() except +

        void setup(int nx, int ny, int nz, double dx, double dy, double dz,
                   double unit_length, double *coordinates, int *ngbs, 
                   double *spin, double *Ms, double *Ms_inv, 
                   double *energy, double *field, int *pins
                   )

        void add_interaction(void * interaction, int int_id)
        void print_interactions_id()


cdef class PyMicroSim(object):
    """
    Wrapper for the C++ MicroSim simulation class
    """
    cdef MicroSim *thisptr
    # Try cinit:
    def __cinit__(self):

        self.thisptr = new MicroSim()

    def __dealloc__(self):
        del self.thisptr

    def setup(self, nx, ny, nz, dx, dy, dz, unit_length,
              double [:, :] coordinates, int [:, :] neighbours,
              double [:] spin, double [:] Ms, double [:] Ms_inv,
              double [:] energy, double [:] field, int [:] pins
              ):

        return self.thisptr.setup(nx, ny, nz, dx, dy, dz, unit_length,
                                  &coordinates[0, 0], &neighbours[0, 0],
                                  &spin[0], &Ms[0], &Ms_inv[0],
                                  &energy[0], &field[0], &pins[0]
                                  )

    def add_interaction(self, Interaction):
        Interaction.add_interaction_to_sim(self)
        self.thisptr.print_interactions_id()

