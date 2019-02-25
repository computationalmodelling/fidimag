# # These header calls Necessary for the Verlet integrator step (not working now)
# cimport numpy as np  # import special compile-time information about numpy
# # Function Prototype for the C function
# ctypedef double (*cfunction)(double t, double * y)
# cdef object f

cdef extern from "nebm_lib.h":

    void compute_tangents_C(double *tangents,
                            double *y,
                            double *energies,
                            int n_dofs_image,
                            int n_images
                            )

    void compute_spring_force_C(double *spring_force,
                                double *y,
                                double *tangents,
                                double *k,
                                int n_images,
                                int n_dofs_image,
                                double * distances
                                )

    void compute_effective_force_C(double * G,
                                   double * tangents,
                                   double * gradientE,
                                   double * spring_force,
                                   int * climbing_image,
                                   int n_images,
                                   int n_dofs_image)

    void normalise(double * a, int n)

    void normalise_images_C(double * y, int n_images, int n_dofs_image)

    void normalise_spins_C(double * y, int n_images, int n_dofs_image)

    double compute_distance_cartesian(double * A, double * B, int n_dofs_image,
                                      int * material, int n_dofs_image_material
                                      )

    void compute_dYdt_C(double * y, double * G, double * dYdt, int * pins,
                        int n_images, int n_dofs_image)

    void compute_dYdt_nc_C(double * y, double * G, double * dYdt, int * pins,
                           int n_images, int n_dofs_image)

    void project_images_C(double * vector, double * y,
                          int n_images, int n_dofs_image,
                          )

    void project_vector_C(double * vector, double * y,
                          int n_dofs_image,
                          )

cdef extern from "nebm_integrators.h":

    double step_Verlet_C(double * forces,
                         double * forces_prev,
                         double * velocities,
                         double * velocities_new,
                         double * y,
                         double t,
                         double h,
                         double mass,
                         int n_images,
                         int n_dofs_image,
                         double (* update_field) (double, double *)
                         )

# Not working
# cdef double cfunction_cb(double t, double [:] y):
#     global f
#     result = f(t, y)
#     return result

# Not working
# def step_Verlet(double [:] forces,
#                 double [:] forces_prev,
#                 double [:] velocities,
#                 double [:] velocities_new,
#                 double [:] y,
#                 double t,
#                 double h,
#                 double mass,
#                 int n_images,
#                 int n_dofs_image,
#                 pythonf
#                 ):
#     global f
#     f = pythonf
# 
#     step_Verlet_C(&forces[0],
#                   &forces_prev[0],
#                   &velocities[0],
#                   &velocities_new[0],
#                   &y[0],
#                   t,
#                   h,
#                   mass,
#                   n_images,
#                   n_dofs_image,
#                   <cfunction> cfunction_cb
#                   )


def compute_tangents(double [:] tangents,
                     double [:] y,
                     double [:] energies,
                     n_dofs_image,
                     n_images
                     ):

    compute_tangents_C(&tangents[0], &y[0], &energies[0],
                       n_dofs_image, n_images
                       )

def compute_spring_force(double [:] spring_force,
                         double [:] y,
                         double [:] tangents,
                         double [:] k,
                         n_images,
                         n_dofs_image,
                         double [:] distances
                         ):

    compute_spring_force_C(&spring_force[0], &y[0], &tangents[0], &k[0],
                           n_images, n_dofs_image,
                           &distances[0],
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

def normalise_clib(double [:] a, n):

    normalise(&a[0], n)

def project_images(double [:]  vector,
                   double [:]  y,
                   n_images, n_dofs_image
                   ):

    project_images_C(&vector[0], &y[0],
                     n_images, n_dofs_image
                     )

def project_vector(double [:]  vector,
                   double [:]  y,
                   n_dofs_image
                   ):

    project_vector_C(&vector[0], &y[0],
                     n_dofs_image
                     )

def normalise_images(double [:] y,
                     n_images, n_dofs_image
                     ):

    normalise_images_C(&y[0], n_images, n_dofs_image)

def normalise_spins(double [:] y,
                    n_images, n_dofs_image
                    ):

    normalise_spins_C(&y[0], n_images, n_dofs_image)

def compute_dYdt(double [:]  y,
                 double [:]  G,
                 double [:]  dYdt,
                 int [:] pins,
                 n_images,
                 n_dofs_image
                 ):

    compute_dYdt_C(&y[0], &G[0], &dYdt[0], &pins[0],
                   n_images, n_dofs_image)

def compute_dYdt_nc(double [:]  y,
                    double [:]  G,
                    double [:]  dYdt,
                    int [:] pins,
                    n_images,
                    n_dofs_image
                    ):

    compute_dYdt_nc_C(&y[0], &G[0], &dYdt[0], &pins[0],
                      n_images, n_dofs_image)
