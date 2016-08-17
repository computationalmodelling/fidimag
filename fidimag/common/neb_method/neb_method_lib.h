#include "math.h"
#define WIDE_PI 3.1415926535897932384626433832795L

void compute_tangents_C(double *ys, double *energy,
                        double *tangents, int image_num, int nodes);

void compute_spring_force_C(double *spring_force,
                            double *y,
                            double *tangents,
                            double k,
                            int n_dofs_image
                            );
