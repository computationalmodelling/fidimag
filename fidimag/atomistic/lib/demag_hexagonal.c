
#include "clib.h"
#include "math.h"
#include "stdlib.h"

void dmi_field_bulk(double *spin, double *field, *coords,
                    double *energy, int n) {

    /* Full calculation of Demag */

	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
        double* rij = malloc(3 * sizeof(double));
        double mrij = 0;

        for (int j = 0; j < n; j++) {

            if(j != i){
                
                for(int k = 0; k < 3; k++) {
                    rij[k] = coords[3 * j + k] - coords[3 * i + k];
                }

                rij_mag = sqrt(rij[0] ** 2 + rij[1] ** 2 + rij[2] ** 2)

                for(int k = 0; k < 3; k++) rij_n[k] = rij[k] / rij_mag;

                mrij = spin[3 * j] * rij[0] + spin[3 * j + 1] * rij[1]
                       + spin[3 * j + 2] * rij[2] ;

                for(int k = 0; k < 3; k++){
                    field[3 * i + k] += (3 * rij_n[k] * mrij - spin(3 * j + k))
                                         / (rij_mag ** 3) ;
                }
            }
        }
        free(rij);
    }
}
