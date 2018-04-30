#include "math.h"
#define WIDE_PI 3.1415926535897932384626433832795L

void normalise_spherical(double *restrict a, int n);

void normalise_images_spherical_C(double *restrict y, int n_images,
                                  int n_dofs_image);

double compute_distance_spherical(double *restrict A, double * B, int n,
                                  int *restrict material, int n_dofs_image_material
                                  );
