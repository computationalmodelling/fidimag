import numpy
cimport numpy as np
np.import_array()

cdef extern from "micro_clib.h":
    void compute_exch_field_micro(double *m, double *field, double *energy, double *Ms_inv,
                          double A, double dx, double dy, double dz,
                          int nx, int ny, int nz, int xperiodic, int yperiodic)
    void dmi_field_bulk(double *m, double *field, double *energy, double *Ms_inv,
                              double D, double dx, double dy, double dz,
                              int nx, int ny, int nz, int xperiodic, int yperiodic)

    void dmi_field_interfacial(double *m, double *field, double *energy, double *Ms_inv,
                    double D, double dx, double dy, double dz,
                    int nx, int ny, int nz, int xperiodic, int yperiodic)



def compute_exchange_field_micro(np.ndarray[double, ndim=1, mode="c"] m,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            np.ndarray[double, ndim=1, mode="c"] energy,
                            np.ndarray[double, ndim=1, mode="c"] Ms_inv,
                            A, dx, dy, dz, nx, ny, nz,
                            xperiodic,
                            yperiodic):

    compute_exch_field_micro(&m[0], &field[0], &energy[0], &Ms_inv[0], A, dx, dy, dz, nx, ny, nz, xperiodic, yperiodic)
    

def compute_dmi_field_bulk(np.ndarray[double, ndim=1, mode="c"] m,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            np.ndarray[double, ndim=1, mode="c"] energy,
                            np.ndarray[double, ndim=1, mode="c"] Ms_inv,
                            D, dx, dy, dz, nx, ny, nz,
                            xperiodic,
                            yperiodic):

    dmi_field_bulk(&m[0], &field[0], &energy[0], &Ms_inv[0], D, dx, dy, dz, nx, ny, nz, xperiodic, yperiodic)
    
def compute_dmi_field_interfacial(np.ndarray[double, ndim=1, mode="c"] m,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            np.ndarray[double, ndim=1, mode="c"] energy,
                            np.ndarray[double, ndim=1, mode="c"] Ms_inv,
                            D, dx, dy, dz, nx, ny, nz,
                            xperiodic,
                            yperiodic):

    dmi_field_interfacial(&m[0], &field[0], &energy[0], &Ms_inv[0], D, dx, dy, dz, nx, ny, nz, xperiodic, yperiodic)