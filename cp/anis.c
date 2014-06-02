#include "clib.h"

void compute_anis(double *spin, double *field, double *Ku, int nxyz) {

	int i;

    #pragma omp parallel for private(i)
	for (i = 0; i < 3*nxyz; i++) {

		field[i] = 2 * Ku[i] * spin[i];

	}

}

double compute_anis_energy(double *spin, double *Ku, int nxyz) {

	int i;

	double energy = 0;

	#pragma omp parallel for private(i) reduction(+:energy)
	for (i = 0; i < 3*nxyz; i++) {

		energy += Ku[i] * spin[i] * spin[i];

	}

	energy = -energy;

	return energy;

}
