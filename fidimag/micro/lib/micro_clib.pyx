import numpy as np
cimport numpy as np

cdef extern from "micro_clib.h":
    void compute_exch_field_micro(double *m, double *field,
                                  double *energy, double *Ms_inv,
                                  double A, double dx, double dy, double dz,
                                  int n, int *ngbs)
    void compute_exch_field_rkky_micro(double *m, double *field, double *energy,
                                  double *Ms_inv, double sigma, int nx, double ny,
                                  double nz, int z_bottom, int z_top)

    void dmi_field(double *m, double *field, double *energy, double *Ms_inv,
                   double *D, int n_dmis,
                   double *dmi_vector,
                   double dx, double dy, double dz,
                   int n, int *ngbs)

    void compute_uniaxial_anis(double *m, double *field,
                               double *energy, double *Ms_inv,
                               double *Ku, double *axis,
                               int nx, int ny, int nz)


    void compute_uniaxial4_anis(double *m, double *field,
                               double *energy, double *Ms_inv,
                               double *K1, double *K2,
                               double *axis,
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


def compute_exchange_field_micro_rkky(double [:] m,
                                      double [:] field,
                                      double [:] energy,
                                      double [:] Ms_inv,
                                 sigma, nx, ny, nz, z_bottom, z_top):

    compute_exch_field_rkky_micro(&m[0], &field[0], &energy[0], &Ms_inv[0], sigma,
                             nx, ny, nz, z_bottom, z_top)


def compute_dmi_field(double [:] m,
                      double [:] field,
                      double [:] energy,
                      double [:] Ms_inv,
                      double [:] D,
                      n_dmis,
                      double [:] dmi_vector,
                      dx, dy, dz,
                      n,
                      int [:, :] ngbs
                      ):

    dmi_field(&m[0], &field[0], &energy[0], &Ms_inv[0],
              &D[0], n_dmis, &dmi_vector[0],
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


def compute_anisotropy4_micro(double [:] m,
                             double [:] field,
                             double [:] energy,
                             double [:] Ms_inv,
                             double [:] K1,
                             double [:] K2,
                             double [:] axis,
                             nx, ny, nz):

    compute_uniaxial4_anis(&m[0], &field[0], &energy[0], &Ms_inv[0],
                          &K1[0], &K2[0], &axis[0], nx, ny, nz)


def compute_skyrmion_number(double [:] m,
                            double [:] charge,
                            nx, ny, nz,
                            int [:, :] ngbs):

    return skyrmion_number(&m[0], &charge[0], nx, ny, nz, &ngbs[0, 0])
