#include "clib.h"


void compute_anis(double *spin, double *field, double *energy,
	double *Ku, double *axis, int nx, int ny, int nz) {
	
	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;

    	#pragma omp parallel for
	for (int i = 0; i < n1; i++) {
	    int j = i+n1;
            int k = i+n2;

           double m_u = spin[i]*axis[i] + spin[j]*axis[j] + spin[k]*axis[k];

		field[i] = 2*Ku[i]*m_u*axis[i];
		field[j] = 2*Ku[i]*m_u*axis[j];
		field[k] = 2*Ku[i]*m_u*axis[k];

		energy[i] = -Ku[i]*(m_u*m_u);

	}

}
