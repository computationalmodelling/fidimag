import numpy
cimport numpy as np
np.import_array()

     

cdef extern from "neb.h":
    void compute_tangents_c(double *ys, double *energy, double *tangents, int image_num, int nodes)
    void compute_dm_dt_c(double *ys, double *heff, double *dm_dt, int *pins, int image_num, int nodes) 


def compute_tangents(np.ndarray[double, ndim=2, mode="c"] ys,
                            np.ndarray[double, ndim=1, mode="c"] energy,
                            np.ndarray[double, ndim=2, mode="c"] tangents,
                            total_image_num,
                            nodes):

    compute_tangents_c(&ys[0,0], &energy[0], &tangents[0,0], total_image_num, nodes)

def compute_dm_dt(np.ndarray[double, ndim=2, mode="c"] ys,
                            np.ndarray[double, ndim=2, mode="c"] heff,
                            np.ndarray[double, ndim=2, mode="c"] dmdt,
                            np.ndarray[int, ndim=1, mode="c"] pins,
                            total_image_num,
                            nodes):

    compute_dm_dt_c(&ys[0,0], &heff[0,0], &dmdt[0,0], &pins[0], total_image_num, nodes)