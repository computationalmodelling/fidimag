#include "math.h"
#define WIDE_PI 3.1415926535897932384626433832795L

void cross_product(double *output, double *A, double *B);

double dot_product(double *A, double *B, int n);

double compute_norm(double *restrict a, int n);

void normalise(double *restrict a, int n);

void compute_tangents_C(double *restrict ys, double *restrict energy,
                        double *restrict tangents, int image_num, int nodes);

void compute_spring_force_C(
    double *restrict spring_force,
    double *restrict y,
    double *restrict tangents,
    double *restrict k,
    int n_images,
    int n_dofs_image,
    double *restrict distances);

void compute_effective_force_C(double *restrict G,
                               double *restrict tangents,
                               double *restrict gradientE,
                               double *restrict spring_force,
                               int *restrict climbing_image,
                               int n_images,
                               int n_dofs_image);

void compute_image_distances(double *restrict distances,
                             double *restrict path_distances,
                             double *restrict y,
                             int n_images,
                             int n_dofs_image,
                             double (*compute_distance)(double *,
                                                        double *,
                                                        int,
                                                        int *,
                                                        int),
                             int *restrict material,
                             int n_dofs_image_material);

void normalise_images_C(double *restrict y, int n_images, int n_dofs_image);

void normalise_spins_C(double *restrict y, int n_images, int n_dofs_image);

void project_images_C(double *restrict vector, double *restrict y,
                      int n_images, int n_dofs_image);

void project_vector_C(double *restrict vector, double *restrict y,
                      int n_dofs_image);

double compute_distance_cartesian(double *restrict A, double *restrict B, int n_dofs_image,
                                  int *restrict material, int n_dofs_image_material);

void compute_dYdt_C(double *restrict y, double *restrict G, double *restrict dYdt,
                    int *restrict pins,
                    int n_images, int n_dofs_image);

void compute_dYdt_nc_C(double *restrict y, double *restrict G, double *restrict dYdt,
                       int *restrict pins,
                       int n_images, int n_dofs_image);
