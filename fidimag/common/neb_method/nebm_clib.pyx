cdef extern from "nebm_lib.h":

    void compute_tangents_C(double *tangents,
                            double *y,
                            double *energies,
                            int n_dofs_image,
                            int n_images
                            )

    void compute_effective_force_C(double * G,
                                   double * tangents,
                                   double * gradientE,
                                   double * spring_force,
                                   int * climbing_image,
                                   int n_images,
                                   int n_dofs_image)

def compute_tangents(double [:] tangents,
                     double [:] y,
                     double [:] energies,
                     n_dofs_image,
                     n_images
                     ):

    compute_tangents_C(&tangents[0], &y[0], &energies[0],
                       n_dofs_image, n_images
                       )

def compute_effective_force(double [:] G,
                            double [:] tangents,
                            double [:] gradientE,
                            double [:] spring_force,
                            int [:] climbing_image,
                            n_images,
                            n_dofs_image
                            ):

    compute_effective_force_C(&G[0], &tangents[0],
                              &gradientE[0], &spring_force[0],
                              &climbing_image[0],
                              n_images, n_dofs_image
                              )
