import numpy
cimport numpy as np
np.import_array()


cdef extern from "nebm_spherical_lib.h":

    void normalise(double * a, int n)

    double compute_distance_spherical(double * A, double * B, int n);
    

cdef extern from "nebm_lib.h":

    void compute_tangents_C(double *tangents,
                            double *y,
                            double *energies,
                            int n_dofs_image,
                            int n_images,
                            void (* normalise)(double *, int)
                            )

    void compute_spring_force_C(double *spring_force,
                                double *y,
                                double *tangents,
                                double k,
                                int n_images,
                                int n_dofs_image,
                                double (* compute_distance)(double *, double *, int)
                                )

    void compute_effective_force_C(double * G,
                                   double * tangents,
                                   double * gradientE,
                                   double * spring_force,
                                   int n_images,
                                   int n_dofs_image)

def compute_tangents(np.ndarray[double, ndim=1, mode="c"] tangents,
                     np.ndarray[double, ndim=1, mode="c"] y,
                     np.ndarray[double, ndim=1, mode="c"] energies,
                     n_dofs_image,
                     n_images
                     ):

    compute_tangents_C(&tangents[0], &y[0], &energies[0],
                       n_dofs_image, n_images,
                       normalise
                       )

def compute_spring_force(np.ndarray[double, ndim=1, mode="c"] spring_force,
                         np.ndarray[double, ndim=1, mode="c"] y,
                         np.ndarray[double, ndim=1, mode="c"] tangents,
                         k,
                         n_images,
                         n_dofs_image,
                         ):

    compute_spring_force_C(&spring_force[0], &y[0], &tangents[0],
                           k, n_images, n_dofs_image,
                           compute_distance_spherical
                           )

def compute_effective_force(np.ndarray[double, ndim=1, mode="c"] G,
                            np.ndarray[double, ndim=1, mode="c"] tangents,
                            np.ndarray[double, ndim=1, mode="c"] gradientE,
                            np.ndarray[double, ndim=1, mode="c"] spring_force,
                            n_images,
                            n_dofs_image
                            ):

    compute_effective_force_C(&G[0], &tangents[0],
                              &gradientE[0], &spring_force[0],
                              n_images, n_dofs_image
                              )
