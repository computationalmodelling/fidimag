#include<math.h>
#define WIDE_PI 3.1415926535897932384626433832795L

void compute_tangents_c(double *ys, double *energy,
                        double *tangents, int image_num, int nodes);

void compute_dm_dt_c(double *ys, double *heff,
                     double *dm_dt, int *pins, int image_num, int nodes) ;
