#include "math.h"
#define WIDE_PI 3.1415926535897932384626433832795L


void cross_product(double * output, double * A, double * B);

double dot_product(double * A, double * B, int n);

double compute_norm(double *a, int n, int scale);

void normalise(double *a, int n);

void compute_tangents_C(double *ys, double *energy,
                        double *tangents, int image_num, int nodes
                        );

void compute_spring_force_C(
        double *spring_force,
        double *y,
        double *tangents,
        double k,
        int n_images,
        int n_dofs_image,
        double (* compute_distance)(double *, double *, int, int *, int),
        int * material,
        int n_dofs_image_material
        );

void compute_effective_force_C(double * G,
                               double * tangents,
                               double * gradientE,
                               double * spring_force,
                               int n_images,
                               int n_dofs_image
                               );

void project_images_C(double * vector, double * y,
                      int n_images, int n_dofs_image
                      );

void normalise_images_C(double * y, int n_images, int n_dofs_image);
