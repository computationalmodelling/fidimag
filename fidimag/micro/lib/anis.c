#include "micro_clib.h"


void compute_uniaxial_anis(double *m, double *field, double *energy, double *Ms_inv, 
	double *Ku, double *axis, int nx, int ny, int nz) {
	
	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;

    #pragma omp parallel for
	for (int i = 0; i < n1; i++) {
		int j = i+n1;
		int k = i+n2;

		if (Ms_inv[i] == 0.0){
			field[i] = 0;
            field[j] = 0;
            field[k] = 0;
            energy[i] = 0;
            continue;
        }

        double m_u = m[i]*axis[i] + m[j]*axis[j] + m[k]*axis[k];

		field[i] = 2*Ku[i]*m_u*Ms_inv[i]*MU0_INV*axis[i];
		field[j] = 2*Ku[i]*m_u*Ms_inv[i]*MU0_INV*axis[j];
		field[k] = 2*Ku[i]*m_u*Ms_inv[i]*MU0_INV*axis[k];

		energy[i] = Ku[i]*(1-m_u*m_u);

	}

}

