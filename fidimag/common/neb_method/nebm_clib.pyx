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

    void normalise(double * a, int n)

    void normalise_images_C(double * y, int n_images, int n_dofs_image)

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
