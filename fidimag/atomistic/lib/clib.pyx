import numpy
cimport numpy as np
np.import_array()


cdef extern from "clib.h":
    void initial_random(int seed)
    void gauss_random_vec(double *x, int n)

    double skyrmion_number(double *spin, double *charge,
                           int nx, int ny, int nz, int *ngbs)

    void compute_guiding_center(double *spin, int nx, int ny, int nz,
                                double *res)

    void compute_px_py_c(double *spin, int nx, int ny, int nz,
                         double *px, double *py)

    void compute_exch_field(double *spin, double *field, double *energy,
                            double Jx, double Jy, double Jz,
                            int *ngbs, int n)
    void compute_exch_field_spatial(double *spin, double *field, double *energy,
                            double *J, int *ngbs, int n)

    double compute_exch_energy(double *spin, double Jx, double Jy, double Jz,
                               int nx, int ny, int nz,
                               int xperiodic, int yperiodic)

    void dmi_field_bulk(double *spin, double *field, double *energy,
                        double *D, int *ngbs, int n)

    void dmi_field_interfacial_atomistic(double *spin, double *field,
                                         double *energy, double D, int *ngbs,
                                         int n, int nneighbours,
                                         double *DMI_vec)

    void demag_full(double *spin, double *field, double *coords,
                    double *energy,
                    double *mu_s, double *mu_s_scale,
                    int n)

    double dmi_energy(double *spin, double D, int nx, int ny, int nz,
                      int xperiodic, int yperiodic)

    void compute_anis(double *spin, double *field, double *energy,
                      double *Ku, double *axis, int n)

    void llg_rhs(double * dm_dt, double * spin,
                 double *h, double *alpha, int *pins,
                 double gamma, int n, int do_precession, double default_c)

    void llg_s_rhs(double * dm_dt, double * spin, double * h,
                   double *alpha, double *chi, double gamma, int n)

    void llg_rhs_jtimes(double *jtn, double *m, double *h,
                        double *mp, double *hp, double *alpha, int *pins,
                        double gamma, int n,
                        int do_precession, double default_c)

    void normalise(double *m, int *pins, int n)

    void compute_stt_field_c(double *spin, double *field,
                             double *jx, double *jy,
                             double dx, double dy, int *ngbs, int n)

    void llg_stt_rhs(double *dm_dt, double *m, double *h, double *h_stt,
                     double *alpha,double beta, double u0, double gamma, int n)

    void llg_stt_cpp(double *dm_dt, double *m, double *h, double *p,
			      double *alpha, int *pins, double *a_J, double beta, double gamma, int n)

    # used for sllg
    void llg_rhs_dw_c(double *m, double *h, double *dm, double *T, double *alpha, double *mu_s_inv, int *pins, double *eta, int n, double gamma, double dt)


def compute_skyrmion_number(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] charge,
                            nx, ny, nz,
                            np.ndarray[int, ndim=2, mode="c"] ngbs
                            ):

    return skyrmion_number(&spin[0], &charge[0], nx, ny, nz, &ngbs[0,0])

def compute_RxRy(np.ndarray[double, ndim=1, mode="c"] spin,
                            nx, ny, nz):

    res = numpy.array([0.0,0.0])

    cdef np.ndarray[double, ndim=1, mode="c"] R = res

    compute_guiding_center(&spin[0], nx, ny, nz, &R[0])

    return res[0], res[1]

def compute_px_py(np.ndarray[double, ndim=1, mode="c"] spin,
                  nx,ny,nz,
                  np.ndarray[double, ndim=1, mode="c"] px,
                  np.ndarray[double, ndim=1, mode="c"] py):

    compute_px_py_c(&spin[0], nx, ny, nz, &px[0], &py[0])

def compute_exchange_field(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            np.ndarray[double, ndim=1, mode="c"] energy,
                            Jx, Jy, Jz,
                            np.ndarray[int, ndim=2, mode="c"] ngbs,
                            n
                            ):

    compute_exch_field(&spin[0], &field[0], &energy[0],
		       Jx, Jy, Jz,
                       &ngbs[0, 0], n)

def compute_exchange_field_spatial(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            np.ndarray[double, ndim=1, mode="c"] energy,
                            np.ndarray[double, ndim=2, mode="c"] J,
                            np.ndarray[int, ndim=2, mode="c"] ngbs,
                            n):

    compute_exch_field_spatial(&spin[0], &field[0], &energy[0],&J[0,0],&ngbs[0, 0], n)


def compute_exchange_energy(np.ndarray[double, ndim=1, mode="c"] spin,
                            Jx, Jy, Jz, nx, ny, nz, xperiodic,yperiodic):

    return compute_exch_energy(&spin[0], Jx, Jy, Jz,
                               nx, ny, nz,
                               xperiodic, yperiodic)


def compute_anisotropy(np.ndarray[double, ndim=1, mode="c"] spin,
                        np.ndarray[double, ndim=1, mode="c"] field,
                        np.ndarray[double, ndim=1, mode="c"] energy,
                        np.ndarray[double, ndim=1, mode="c"] Ku,
                        np.ndarray[double, ndim=1, mode="c"] axis,
                        n):
    compute_anis(&spin[0], &field[0], &energy[0], &Ku[0],
                 &axis[0], n)


def compute_dmi_field(np.ndarray[double, ndim=1, mode="c"] spin,
                      np.ndarray[double, ndim=1, mode="c"] field,
                      np.ndarray[double, ndim=1, mode="c"] energy,
                      np.ndarray[double, ndim=2, mode="c"] D,
                      np.ndarray[int, ndim=2, mode="c"] ngbs,
                      n):
    dmi_field_bulk(&spin[0], &field[0], &energy[0], &D[0,0], &ngbs[0, 0], n)


def compute_dmi_field_interfacial(np.ndarray[double, ndim=1, mode="c"] spin,
                                  np.ndarray[double, ndim=1, mode="c"] field,
                                  np.ndarray[double, ndim=1, mode="c"] energy,
                                  D,
                                  np.ndarray[int, ndim=2, mode="c"] ngbs,
                                  n, nneighbours,
                                  np.ndarray[double, ndim=1, mode="c"] DMI_vec,
                                  ):
    dmi_field_interfacial_atomistic(&spin[0], &field[0], &energy[0],
                                    D, &ngbs[0, 0], n,
                                    nneighbours,
                                    &DMI_vec[0]
                                    )

def compute_demag_full(np.ndarray[double, ndim=1, mode="c"] spin,
                       np.ndarray[double, ndim=1, mode="c"] field,
                       np.ndarray[double, ndim=1, mode="c"] energy,
                       np.ndarray[double, ndim=2, mode="c"] coords,
                       np.ndarray[double, ndim=1, mode="c"] mu_s,
                       np.ndarray[double, ndim=1, mode="c"] mu_s_scale,
                       n
                       ):
    demag_full(&spin[0], &field[0], &energy[0],
               &coords[0, 0], &mu_s[0], &mu_s_scale[0], n)

def compute_dmi_energy(np.ndarray[double, ndim=1, mode="c"] spin,
                        D, nx, ny, nz,
                        xperiodic,yperiodic):
    return dmi_energy(&spin[0], D, nx, ny, nz, xperiodic,yperiodic)



def compute_llg_rhs(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] alpha,
                np.ndarray[int, ndim=1, mode="c"] pins,
                gamma, n, do_precession, default_c):
    llg_rhs(&dm_dt[0], &spin[0], &field[0], &alpha[0], &pins[0],
            gamma, n, do_precession, default_c)


def compute_llg_jtimes(np.ndarray[double, ndim=1, mode="c"] jtn,
                np.ndarray[double, ndim=1, mode="c"] m,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] mp,
                np.ndarray[double, ndim=1, mode="c"] field_p,
                np.ndarray[double, ndim=1, mode="c"] alpha,
                np.ndarray[int, ndim=1, mode="c"] pins,
                gamma, n, do_precession, default_c):
    llg_rhs_jtimes(&jtn[0], &m[0], &field[0], &mp[0], &field_p[0],
                   &alpha[0], &pins[0], gamma, n, do_precession, default_c)

def compute_llg_s_rhs(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] alpha,
                np.ndarray[double, ndim=1, mode="c"] chi,
                gamma, n):
    llg_s_rhs(&dm_dt[0], &spin[0], &field[0], &alpha[0], &chi[0], gamma, n)


def compute_stt_field(np.ndarray[double, ndim=1, mode="c"] spin,
                      np.ndarray[double, ndim=1, mode="c"] field,
                      np.ndarray[double, ndim=1, mode="c"] jx,
                      np.ndarray[double, ndim=1, mode="c"] jy,
                      dx, dy,
                      np.ndarray[int, ndim=2, mode="c"] ngbs,
                      n
                      ):
    compute_stt_field_c(&spin[0], &field[0], &jx[0], &jy[0],
                        dx, dy, &ngbs[0, 0], n)

def compute_llg_stt_rhs(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] field_stt,
                np.ndarray[double, ndim=1, mode="c"] alpha,
                beta, u0, gamma, n):
    llg_stt_rhs(&dm_dt[0], &spin[0], &field[0], &field_stt[0],
                &alpha[0], beta, u0, gamma, n)



def compute_llg_stt_cpp(np.ndarray[double, ndim=1, mode="c"] dm_dt,
                np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] p,
                np.ndarray[double, ndim=1, mode="c"] alpha,
		            np.ndarray[int, ndim=1, mode="c"] pin,
                np.ndarray[double, ndim=1, mode="c"] a_J,
                beta, gamma, n):
    llg_stt_cpp(&dm_dt[0], &spin[0], &field[0], &p[0], 
                &alpha[0], &pin[0], &a_J[0], beta, gamma, n)


def normalise_spin(np.ndarray[double, ndim=1, mode="c"] spin,
              np.ndarray[int, ndim=1, mode="c"] pins, n):
    normalise(&spin[0], &pins[0], n)


def init_random(seed):
    initial_random(seed);

def random_number_array(np.ndarray[double, ndim=1, mode="c"] v):
    
    cdef int n = len(v)

    gauss_random_vec(&v[0], n)


def compute_llg_rhs_dw(np.ndarray[double, ndim=1, mode="c"] dm,
                np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] T,
                np.ndarray[double, ndim=1, mode="c"] alpha,
                np.ndarray[double, ndim=1, mode="c"] mu_s_inv,
                np.ndarray[double, ndim=1, mode="c"] eta,
                np.ndarray[int, ndim=1, mode="c"] pin, n, gamma, dt):
    llg_rhs_dw_c(&spin[0], &field[0], &dm[0], &T[0], &alpha[0],  &mu_s_inv[0], &pin[0], &eta[0], n, gamma, dt)
