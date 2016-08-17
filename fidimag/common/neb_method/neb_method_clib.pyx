import numpy
cimport numpy as np
np.import_array()


cdef extern from "neb_method_lib.h":
    void compute_tangents_C(double *tangents, 
                            double *y,
                            double *energies,
                            int n_dofs_image,
                            int n_images
                            )

    void compute_spring_force_C(double *spring_force,
                                double *y,
                                double *tangents,
                                double k,
                                int n_dofs_image
                                )


def compute_tangents(np.ndarray[double, ndim=1, mode="c"] tangents,
                     np.ndarray[double, ndim=1, mode="c"] y,
                     np.ndarray[double, ndim=1, mode="c"] energies,
                     n_dofs_image,
                     n_images
                     ):

    compute_tangents_C(&tangents[0], &y[0], &energies[0],
                       n_dofs_image, n_images
                       )

def compute_spring_force(np.ndarray[double, ndim=1, mode="c"] spring_force,
                         np.ndarray[double, ndim=1, mode="c"] y,
                         np.ndarray[double, ndim=1, mode="c"] tangents,
                         k,
                         n_dofs_image,
                         ):

    compute_spring_force_C(&spring_force[0], &y[0], &tangents[0],
                           k, n_dofs_image
                           )
