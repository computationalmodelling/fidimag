#include "math.h"
#define WIDE_PI 3.1415926535897932384626433832795L

double compute_distance_cartesian(double * A, double * B, int n);

void compute_dYdt_C(double * y, double * G, double * dYdt,
                    int * pins,
                    int n_images, int n_dofs_image);
