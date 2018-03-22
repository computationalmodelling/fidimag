import numpy
cimport numpy as np
np.import_array()

cdef extern from "fidimag_random.h":
    ctypedef struct mt19937_state:
        pass

    mt19937_state *create_mt19937_state()
    void finalize_mt19937_state(mt19937_state *state)

    void initial_rng_mt19973(mt19937_state *state, int seed)
    double random_double_half_open(mt19937_state *state)
    void gauss_random_vector(mt19937_state *state, double *x, int n)
    void uniform_random_sphere(mt19937_state *state, double *spin, int n)

cdef extern from "time.h":
    ctypedef int time_t
    time_t time(time_t *timer)

cdef extern from "clib.h":
    void run_step_mc(mt19937_state *state, double *spin, double *new_spin,
                     int *ngbs, int *nngbs, int n_ngbs,
                     double J, double J1, double D, double D1,
                     double *h, double Kc, int n, double T, int hexagonal_mesh)

    # -------------------------------------------------------------------------

    double skyrmion_number(double *spin, double *charge,
                           int nx, int ny, int nz, int *ngbs, int n_ngbs)


    double skyrmion_number_BergLuscher(double *spin, double *charge,
                                       int nx, int ny, int nz, int *ngbs, int n_ngbs)

    void compute_guiding_center(double *spin, int nx, int ny, int nz, int nx_start,
                                int nx_stop, int ny_start, int ny_stop,
                                double *res)

    void compute_px_py_c(double *spin, int nx, int ny, int nz,
                         double *px, double *py)

    # -------------------------------------------------------------------------

    void compute_exch_field(double *spin, double *field, double *energy,
                            double Jx, double Jy, double Jz,
                            int *ngbs, int n, int n_ngbs)
    void compute_exch_field_spatial(double *spin, double *field, double *energy,
                            double *J, int *ngbs, int n, int n_ngbs)

    double compute_exch_energy(double *spin, double Jx, double Jy, double Jz,
                               int nx, int ny, int nz,
                               int xperiodic, int yperiodic)

    void compute_full_exch_field(double *spin, double *field, double *energy,
					      	     double *J, int *ngbs, int n, int n_ngbs,
                                 int n_shells, int *n_ngbs_shell,
                                 int *sum_ngbs_shell
                                 )

    # -------------------------------------------------------------------------

    void dmi_field_bulk(double *spin, double *field, double *energy,
                        double *D, int *ngbs, int n, int n_ngbs)

    void dmi_field_interfacial_atomistic(double *spin, double *field,
                                         double *energy, double D, int *ngbs,
                                         int n, int n_ngbs, int n_ngbs_dmi,
                                         double *DMI_vec)

    double dmi_energy(double *spin, double D, int nx, int ny, int nz,
                      int xperiodic, int yperiodic)

    # -------------------------------------------------------------------------

    void demag_full(double *spin, double *field, double *coords,
                    double *energy,
                    double *mu_s, double *mu_s_scale,
                    int n)

    # -------------------------------------------------------------------------

    void compute_anis(double *spin, double *field, double *energy,
                      double *Ku, double *axis, int n)

    void compute_anis_cubic(double *spin, double *field, double *energy, double *Kc, int n)

    # -------------------------------------------------------------------------

    void normalise(double *m, int *pins, int n)

    # used for sllg
    void llg_rhs_dw_c(double *m, double *h, double *dm, double *T, double *alpha, double *mu_s_inv, int *pins, double *eta, int n, double gamma, double dt)

    # -------------------------------------------------------------------------


    void sd_update_spin (double *spin, double *spin_last,
                         double *mxH, double *mxmxH, double *mxmxH_last, double *tau,
                         int* pins, int n)

    void sd_compute_step (double *spin, double *spin_last, double *field, double *scale,
                          double *mxH, double *mxmxH, double *mxmxH_last, double *tau,
                          int *pins, int n, int counter, double tmin, double tmax)



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def compute_skyrmion_number(double [:] spin,
                            double [:] charge,
                            nx, ny, nz,
                            int [:, :] ngbs,
                            n_ngbs
                            ):

    return skyrmion_number(&spin[0], &charge[0], nx, ny, nz, &ngbs[0,0], n_ngbs)


def compute_skyrmion_number_BergLuscher(double [:] spin,
                                        double [:] charge,
                                        nx, ny, nz,
                                        int [:, :] ngbs,
                                        n_ngbs
                                        ):

    return skyrmion_number_BergLuscher(&spin[0], &charge[0], nx, ny, nz, &ngbs[0,0], n_ngbs)

def compute_RxRy(double [:] spin,
                 nx, ny, nz, nx_start=0, nx_stop=-1, ny_start=0, ny_stop=-1):

    res = numpy.array([0.0,0.0])
    if nx_stop < 0 or nx_stop > nx:
        nx_stop = nx
    if ny_stop < 0 or ny_stop > ny:
        ny_stop = ny

    cdef double [:] R = res

    compute_guiding_center(&spin[0], nx, ny, nz, nx_start, nx_stop, ny_start, ny_stop, &R[0])

    return res[0], res[1]

def compute_px_py(double [:] spin,
                  nx,ny,nz,
                  double[:] px,
                  double[:] py):

    compute_px_py_c(&spin[0], nx, ny, nz, &px[0], &py[0])

# -------------------------------------------------------------------------

def compute_exchange_field(double [:] spin,
                           double [:] field,
                           double [:] energy,
                           Jx, Jy, Jz,
                           int [:, :] ngbs,
                           n, n_ngbs
                           ):

    compute_exch_field(&spin[0], &field[0], &energy[0], Jx, Jy, Jz,
                       &ngbs[0, 0], n, n_ngbs)

def compute_exchange_field_spatial(double [:] spin,
                                   double [:] field,
                                   double [:] energy,
                                   double [:, :] J,
                                   int [:, :] ngbs,
                                   n, n_ngbs):

    compute_exch_field_spatial(&spin[0], &field[0], &energy[0],&J[0,0],&ngbs[0, 0], n, n_ngbs)


def compute_exchange_energy(double [:] spin,
                            Jx, Jy, Jz, nx, ny, nz, xperiodic,yperiodic):

    return compute_exch_energy(&spin[0], Jx, Jy, Jz,
                               nx, ny, nz,
                               xperiodic, yperiodic)


def compute_full_exchange_field(double [:] spin,
                                double [:] field,
                                double [:] energy,
                                double [:] J,
                                int [:, :] ngbs,
                                n, n_ngbs, n_shells,
                                int [:] n_ngbs_shell,
                                int [:] sum_ngbs_shell
                                ):

    compute_full_exch_field(&spin[0], &field[0], &energy[0], &J[0],
                            &ngbs[0, 0], n, n_ngbs, n_shells,
                            &n_ngbs_shell[0], &sum_ngbs_shell[0])

# -------------------------------------------------------------------------

def compute_anisotropy(double [:] spin, double [:] field, double [:] energy,
                       double [:] Ku, double [:] axis, n):
    compute_anis(&spin[0], &field[0], &energy[0], &Ku[0], &axis[0], n)

def compute_anisotropy_cubic(double [:] spin, double [:] field, double [:] energy,
                             double [:] Kc, n):
    compute_anis_cubic(&spin[0], &field[0], &energy[0], &Kc[0], n)

# -----------------------------------------------------------------------------

def compute_dmi_field(double [:] spin,
                      double [:] field,
                      double [:] energy,
                      double [:, :] D,
                      int [:, :] ngbs,
                      n, n_ngbs):
    dmi_field_bulk(&spin[0], &field[0], &energy[0], &D[0,0], &ngbs[0, 0], n, n_ngbs)


def compute_dmi_field_interfacial(double [:] spin,
                                  double [:] field,
                                  double [:] energy,
                                  D,
                                  int [:, :] ngbs,
                                  n, n_ngbs, n_ngbs_dmi,
                                  double [:] DMI_vec,
                                  ):
    dmi_field_interfacial_atomistic(&spin[0], &field[0], &energy[0],
                                    D, &ngbs[0, 0], n,
                                    n_ngbs, n_ngbs_dmi,
                                    &DMI_vec[0]
                                    )

def compute_dmi_energy(np.ndarray[double, ndim=1, mode="c"] spin,
                        D, nx, ny, nz,
                        xperiodic,yperiodic):
    return dmi_energy(&spin[0], D, nx, ny, nz, xperiodic,yperiodic)

# -------------------------------------------------------------------------

def compute_demag_full(double [:] spin,
                       double [:] field,
                       double [:] energy,
                       double [:, :] coords,
                       double [:] mu_s,
                       double [:] mu_s_scale,
                       n
                       ):
    demag_full(&spin[0], &field[0], &energy[0],
               &coords[0, 0], &mu_s[0], &mu_s_scale[0], n)

# -------------------------------------------------------------------------

def normalise_spin(np.ndarray[double, ndim=1, mode="c"] spin,
              np.ndarray[int, ndim=1, mode="c"] pins, n):
    normalise(&spin[0], &pins[0], n)

def compute_llg_rhs_dw(double [:] dm,
                       double [:] spin,
                       double [:] field,
                       double [:] T,
                       double [:] alpha,
                       double [:] mu_s_inv,
                       double [:] eta,
                       int [:] pin, n, gamma, dt):
    llg_rhs_dw_c(&spin[0], &field[0], &dm[0], &T[0], &alpha[0],  &mu_s_inv[0], &pin[0], &eta[0], n, gamma, dt)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


cdef class rng_mt19937(object):
    cdef mt19937_state *_c_state
    cdef public int seed
    def __init__(self, seed=None):
        if seed:
            self.seed = int(seed)
        else:
            self.seed = time(NULL)

        self._c_state = create_mt19937_state()
        if self._c_state is NULL:
            raise MemoryError()

        initial_rng_mt19973(self._c_state, self.seed)

    def set_seed(self, seed):
        self.seed = int(seed)
        initial_rng_mt19973(self._c_state, self.seed)

    def random(self):
        """
            return a random number in [0,1)
        """
        return random_double_half_open(self._c_state)

    def fill_vector_gaussian(self, np.ndarray[np.float64_t, ndim=1] vector):
        gauss_random_vector(self._c_state, &vector[0], vector.shape[0])

    def fill_vector_uniform_sphere(self, np.ndarray[double, ndim=1, mode="c"] spin, n):
        uniform_random_sphere(self._c_state,&spin[0], n)

    def __dealloc__(self):
        if self._c_state is not NULL:
            finalize_mt19937_state(self._c_state)
            self._c_state = NULL

cdef class monte_carlo(object):
    cdef mt19937_state *_c_state

    def __init__(self, seed=-1):

        self._c_state = create_mt19937_state()
        if self._c_state is NULL:
            raise MemoryError()

        initial_rng_mt19973(self._c_state, seed)

    def set_seed(self, seed):
        initial_rng_mt19973(self._c_state, seed)

    def run_step(self,
                 double [:] spin,
                 double [:] new_spin,
                 int [:, :] ngbs,
                 int [:, :] nngbs,
                 n_ngbs,
                 J, J1, D, D1, np.ndarray[double, ndim=1, mode="c"] h,
                 Kc, n, T, hexagonal_mesh):

        run_step_mc(self._c_state, &spin[0], &new_spin[0],
                    &ngbs[0,0], &nngbs[0,0], n_ngbs,
                    J, J1, D, D1, &h[0], Kc, n, T, hexagonal_mesh)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def compute_sd_spin(double [:] spin,
                    double [:] spin_last,
                    double [:] mxH,
                    double [:] mxmxH,
                    double [:] mxmxH_last,
                    double [:] tau,
                    int [:] pins,
                    n):

    sd_update_spin(&spin[0], &spin_last[0], &mxH[0],
                   &mxmxH[0], &mxmxH_last[0], &tau[0], &pins[0], n
                   )

def compute_sd_step(double [:] spin,
                    double [:] spin_last,
                    double [:] field,
                    double [:] scale,
                    double [:] mxH,
                    double [:] mxmxH,
                    double [:] mxmxH_last,
                    double [:] tau,
                    int [:] pins,
                    n, counter, tmin, tmax):

    sd_compute_step(&spin[0], &spin_last[0], &field[0], &scale[0], &mxH[0],
                    &mxmxH[0], &mxmxH_last[0], &tau[0], &pins[0],
                    n, counter, tmin, tmax
                    )
