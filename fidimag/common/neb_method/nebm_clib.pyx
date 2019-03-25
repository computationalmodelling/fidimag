cimport nebm_clib

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

    void compute_image_distances(double * distances, double * path_distances,
                                 double * y, int n_images, int n_dofs_image,
                                 double (* compute_distance)(double *,
                                                             double *,
                                                             int,
                                                             int *,
                                                             int),
                                 int *  material, int n_dofs_image_material
                                 )


    void compute_effective_force_C(double * G,
                                   double * tangents,
                                   double * gradientE,
                                   double * spring_force,
                                   int * climbing_image,
                                   int n_images,
                                   int n_dofs_image)

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



cdef extern from "nebm_geodesic_lib.h":
    double compute_geodesic_GreatCircle(double * A, double * B,
                                     int n_dofs_image,
                                     int * material,
                                     int n_dofs_image_material
                                     )

    double compute_geodesic_Vincenty(double * A, double * B,
                                     int n_dofs_image,
                                     int * material,
                                     int n_dofs_image_material
                                     )

def geodesic_distance_Vincenty(double [:] A,
                             double [:] B,
                             n_dofs_image,
                             int [:] material,
                             n_dofs_image_material
                             ):

    return compute_geodesic_Vincenty(&A[0], &B[0], n_dofs_image,
                                     &material[0], n_dofs_image_material
                                     )

def geodesic_distance_GreatCircle(double [:] A,
                                double [:] B,
                                n_dofs_image,
                                int [:] material,
                                n_dofs_image_material
                                ):
    return compute_geodesic_GreatCircle(&A[0], &B[0], n_dofs_image,
                                        &material[0], n_dofs_image_material
                                        )

def image_distances_GreatCircle(double [:] distances,
                                double [:] path_distances,
                                double [:] y,
                                int n_images,
                                int n_dofs_image,
                                int [:] material,
                                int n_dofs_image_material
                                ):
    return compute_image_distances(&distances[0],
                                   &path_distances[0],
                                   &y[0],
                                   n_images,
                                   n_dofs_image,
                                   compute_geodesic_GreatCircle,
                                   &material[0],
                                   n_dofs_image_material
                                   )


cdef extern from "nebm_spherical_lib.h":
    void normalise_spherical(double * a, int n)
    void normalise_images_spherical_C(double * y, int n_images,
                                      int n_dofs_image)
    double compute_distance_spherical(double * A, double * B, int n,
                                      int * material, int n_dofs_image_material
                                      )

def normalise_images_spherical(double [:] y,
                     n_images,
                     n_dofs_image
                     ):
    normalise_images_spherical_C(&y[0], n_images, n_dofs_image)

def image_distances_Spherical(double [:] distances,
                              double [:] path_distances,
                              double [:] y,
                              int n_images,
                              int n_dofs_image,
                              int [:] material,
                              int n_dofs_image_material
                              ):
    return compute_image_distances(&distances[0],
                                   &path_distances[0],
                                   &y[0],
                                   n_images,
                                   n_dofs_image,
                                   compute_distance_spherical,
                                   &material[0],
                                   n_dofs_image_material
                                   )
