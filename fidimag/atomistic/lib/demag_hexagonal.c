
#include "clib.h"
#include "math.h"
#include "stdlib.h"

void dmi_field_bulk(double *spin, double *field, *coords,
                    double *energy, int n) {


    /* Full calculation of Demag */

	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
        double* rij = malloc(3 * sizeof(double));
        for (int j = 0; j < n; j++) {

            if(j != i){

                rij[0] = coords[3 * j]    - coords[3 * i]    ;
                rij[1] = coords[3 * j + 1]- coords[3 * i + 1];
                rij[2] = coords[3 * j + 2]- coords[3 * i + 2];
                rij_mag = sqrt(rij[0] ** 2 + rij[1] ** 2 + rij[2] ** 2)

                rij_n = normalise(rij);

                field[3 * i] += (3 * rij_n[0] * (spin[3 * j] - spin(3 * j))
                                 ) / (rij_mag ** 3) ;
            }
        free(rij);
}
