#include "clib.h"

void compute_anis(double *spin, double *field, double Dx, double Dy, double Dz,
		int nxyz) {

	int i, j, k;

	for (i = 0; i < nxyz; i++) {

		j = i + nxyz;
		k = j + nxyz;

		field[i] = 2 * Dx * spin[i];
		field[j] = 2 * Dy * spin[j];
		field[k] = 2 * Dz * spin[k];

	}

}

double compute_anis_energy(double *spin, double Dx, double Dy, double Dz, int nxyz) {

	int i, j, k;

	double energy=0;

	for (i = 0; i < nxyz; i++) {

		j = i + nxyz;
		k = j + nxyz;

		energy += Dx * spin[i]*spin[i];
		energy += Dy * spin[j]*spin[j];
		energy += Dz * spin[k]*spin[k];

	}

	return energy;

}
