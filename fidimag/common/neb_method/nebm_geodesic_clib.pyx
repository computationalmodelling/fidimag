import numpy
cimport numpy as np
np.import_array()


cdef extern from "nebm_geodesic_lib.h":

    double compute_distance_geodesic(double * A, double * B, int n_dofs_image,
                                     int * material, int n_dofs_image_material
                                     )

cdef extern from "nebm_lib.h":

    void compute_spring_force_C(double *spring_force,
                                double *y,
                                double *tangents,
                                double k,
                                int n_images,
                                int n_dofs_image,
                                double (* compute_distance)(double *, double *, int,
                                                            int *, int
                                                            ),
                                int * material, int n_dofs_image_material
                                )


def compute_spring_force(np.ndarray[double, ndim=1, mode="c"] spring_force,
                         np.ndarray[double, ndim=1, mode="c"] y,
                         np.ndarray[double, ndim=1, mode="c"] tangents,
                         k,
                         n_images,
                         n_dofs_image,
                         np.ndarray[int, ndim=1, mode="c"] material,
                         n_dofs_image_material
                         ):

    compute_spring_force_C(&spring_force[0], &y[0], &tangents[0],
                           k, n_images, n_dofs_image,
                           compute_distance_geodesic,
                           &material[0], n_dofs_image_material
                           )

def geodesic_distance(np.ndarray[double, ndim=1, mode="c"] A,
                      np.ndarray[double, ndim=1, mode="c"] B,
                      n_dofs_image,
                      np.ndarray[int, ndim=1, mode="c"] material,
                      n_dofs_image_material
                      ):

    return compute_distance_geodesic(&A[0], &B[0], n_dofs_image,
                                     &material[0], n_dofs_image_material
                                     )

