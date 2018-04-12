import numpy

cdef extern from "micro_clib.h":
    void compute_exch_field_micro(double *m, double *field,
                                  double *energy, double *Ms_inv,
                                  double A, double dx, double dy, double dz,
                                  int n, int *ngbs)

    void dmi_field(double *m, double *field, double *energy, double *Ms_inv,
                   double *D, double dmi_vector[18], int n_dmi_ngbs,
                   double dx, double dy, double dz, 
                   int n, int *ngbs)

    void compute_uniaxial_anis(double *m, double *field,
                               double *energy, double *Ms_inv, 
                               double *Ku, double *axis,
                               int nx, int ny, int nz)
    
    double skyrmion_number(double *m, double *charge,
                           int nx, int ny, int nz, int *ngbs)


def compute_exchange_field_micro(double [:] m,
                                 double [:] field,
                                 double [:] energy,
                                 double [:] Ms_inv,
                                 A, dx, dy, dz, n,
                                 int [:, :] ngbs):

    compute_exch_field_micro(&m[0], &field[0], &energy[0], &Ms_inv[0], A,
                             dx, dy, dz, n, &ngbs[0, 0])
    

def compute_dmi_field(double [:] m,
                      double [:] field,
                      double [:] energy,
                      double [:] Ms_inv,
                      double [:] D,
                      double [:] dmi_vector,
                      n_dmi_ngbs,
                      dx, dy, dz,
                      n,
                      int [:, :] ngbs
                      ):

    dmi_field(&m[0], &field[0], &energy[0], &Ms_inv[0], 
              &D[0], &dmi_vector[0], n_dmi_ngbs,
              dx, dy, dz, n, &ngbs[0, 0])
    

def compute_anisotropy_micro(double [:] m,
                             double [:] field,
                             double [:] energy,
                             double [:] Ms_inv,
                             double [:] Ku,
                             double [:] axis,
                             nx, ny, nz):

    compute_uniaxial_anis(&m[0], &field[0], &energy[0], &Ms_inv[0],
                          &Ku[0], &axis[0], nx, ny, nz)

def compute_skyrmion_number(double [:] m,
                            double [:] charge,
                            nx, ny, nz,
                            int [:, :] ngbs):

    return skyrmion_number(&m[0], &charge[0], nx, ny, nz, &ngbs[0, 0])
