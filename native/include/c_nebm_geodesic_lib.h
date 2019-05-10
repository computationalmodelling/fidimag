#include "math.h"

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
