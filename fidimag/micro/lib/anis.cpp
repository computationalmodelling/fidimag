#include "micro_clib.h"


void compute_uniaxial_anis(double * m, double * field, double * energy, double * Ms_inv, 
	double * Ku, double * axis, int nx, int ny, int nz) {
	
	int n = nx * ny * nz;

    #pragma omp parallel for
	for (int i = 0; i < n; i++) {
		int j = 3 * i;

		if (Ms_inv[i] == 0.0){
			field[j] = 0;
            field[j + 1] = 0;
            field[j + 2] = 0;
            energy[i] = 0;
            continue;
        }

        double m_u = m[j] * axis[j] + m[j + 1] * axis[j + 1] + m[j + 2] * axis[j + 2];

		field[j]     = 2 * Ku[i] * m_u * Ms_inv[i] * MU0_INV * axis[j];
		field[j + 1] = 2 * Ku[i] * m_u * Ms_inv[i] * MU0_INV * axis[j + 1];
		field[j + 2] = 2 * Ku[i] * m_u * Ms_inv[i] * MU0_INV * axis[j + 2];

		energy[i] = Ku[i] * (1 - m_u * m_u);

	}

}

