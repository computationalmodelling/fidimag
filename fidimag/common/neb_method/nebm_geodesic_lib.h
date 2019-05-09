#include "math.h"
#define WIDE_PI 3.1415926535897932384626433832795L

double compute_geodesic_GreatCircle(double * A, double * B,
                                    int n_dofs_image,
                                    int * material,
                                    int n_dofs_image_material
                                    );

double compute_geodesic_Vincenty(double * A, double * B,
                                 int n_dofs_image,
                                 int * material,
                                 int n_dofs_image_material
                                 );
