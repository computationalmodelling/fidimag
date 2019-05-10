#include "math.h"

void compute_tangents_C(double * ys, double * energy,
                        double * tangents, int image_num, int nodes
                        );

void compute_spring_force_C(
        double * spring_force,
        double * y,
        double * tangents,
        double * k,
        int n_images,
        int n_dofs_image,
        double * distances
        );

void compute_effective_force_C(double * G,
                               double * tangents,
                               double * gradientE,
                               double * spring_force,
                               int    * climbing_image,
                               int n_images,
                               int n_dofs_image
                               );

void compute_image_distances(double * distances,
                             double * path_distances,
                             double * y,
                             int n_images,
                             int n_dofs_image,
                             double (* compute_distance)(double *,
                                                         double *,
                                                         int,
                                                         int *,
                                                         int),
                             int *  material,
                             int n_dofs_image_material
                             );

void normalise_images_C(double * y, int n_images, int n_dofs_image);

void normalise_spins_C(double * y, int n_images, int n_dofs_image);

void project_images_C(double * vector, double * y,
                      int n_images, int n_dofs_image
                      );

void project_vector_C(double * vector, double * y,
                      int n_dofs_image
                      );

double compute_distance_cartesian(double * A, double * B, int n_dofs_image,
                                  int * material, int n_dofs_image_material
                                  );

void compute_dYdt_C(double * y, double * G, double * dYdt,
                    int * pins,
                    int n_images, int n_dofs_image);

void compute_dYdt_nc_C(double * y, double * G, double * dYdt,
                       int * pins,
                       int n_images, int n_dofs_image);
