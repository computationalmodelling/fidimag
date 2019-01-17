#include "math.h"
#define WIDE_PI 3.1415926535897932384626433832795L

double compute_distance_cartesian(double *restrict A, double *restrict B, int n_dofs_image,
                                  int *restrict material, int n_dofs_image_material
                                  );

void compute_dYdt_C(double *restrict y, double *restrict G, double *restrict dYdt,
                    int *restrict pins,
                    int n_images, int n_dofs_image);

void compute_dYdt_nc_C(double *restrict y, double *restrict G, double *restrict dYdt,
                       int *restrict pins,
                       int n_images, int n_dofs_image);


void project_images_C(double *restrict vector, double *restrict y,
                      int n_images, int n_dofs_image
                      );
