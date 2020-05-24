#include "micro_clib.h"


void compute_uniaxial_anis(double *restrict m, double *restrict field, double *restrict energy, double *restrict Ms_inv, 
	double *restrict Ku, double *restrict axis, int nx, int ny, int nz) {
	
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


void compute_uniaxial4_anis(double *restrict m, double *restrict field, double *restrict energy, double *restrict Ms_inv, 
    double *restrict K1, double *restrict K2, double *restrict axis, int nx, int ny, int nz) {
    
    int n = nx * ny * nz;

    // Follows calculation of OOMMF extension by Hans and Richard Boardman
    // http://www.soton.ac.uk/~fangohr/software/oxs_uniaxial4/download/uniaxialanisotropy4.cc

    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        int j = 3 * i;

        if (Ms_inv[i] == 0.0) {
            field[j] = 0;
            field[j + 1] = 0;
            field[j + 2] = 0;
            energy[i] = 0;
            continue;
        }

        double k1 = K1[i];
        double k2 = K2[i];

        double field_mult1 = MU0_INV * 2.0 * k1 * Ms_inv[i];
        double field_mult2 = MU0_INV * 4.0 * k2 * Ms_inv[i];
        double m_dot_u = m[j] * axis[j] + m[j + 1] * axis[j + 1] + m[j + 2] * axis[j + 2];

        if (k1 <= 0) {
            field[j + 0] = (field_mult1*m_dot_u) * axis[j + 0] + (field_mult2 * m_dot_u*m_dot_u*m_dot_u) * axis[j + 0];
            field[j + 1] = (field_mult1*m_dot_u) * axis[j + 1] + (field_mult2 * m_dot_u*m_dot_u*m_dot_u) * axis[j + 1];
            field[j + 2] = (field_mult1*m_dot_u) * axis[j + 2] + (field_mult2 * m_dot_u*m_dot_u*m_dot_u) * axis[j + 2];
            energy[i] = -k1*m_dot_u*m_dot_u - k2*m_dot_u*m_dot_u*m_dot_u*m_dot_u;
        }

        else {
            double u_x_m[3];
            u_x_m[0] = cross_x(axis[j], axis[j+1], axis[j+2], m[j], m[j+1], m[j+2]);
            u_x_m[1] = cross_y(axis[j], axis[j+1], axis[j+2], m[j], m[j+1], m[j+2]);
            u_x_m[2] = cross_z(axis[j], axis[j+1], axis[j+2], m[j], m[j+1], m[j+2]);
            double u_x_m_mag2 = u_x_m[1]*u_x_m[1] + u_x_m[1]*u_x_m[1] + u_x_m[2]*u_x_m[2];
            field[j + 0] = (field_mult1*m_dot_u) * axis[j + 0] + (field_mult2*m_dot_u*m_dot_u*m_dot_u) * axis[j + 0];
            field[j + 1] = (field_mult1*m_dot_u) * axis[j + 1] + (field_mult2*m_dot_u*m_dot_u*m_dot_u) * axis[j + 1];
            field[j + 2] = (field_mult1*m_dot_u) * axis[j + 2] + (field_mult2*m_dot_u*m_dot_u*m_dot_u) * axis[j + 2];
            energy[i] = (k1 + 2*k2)*u_x_m_mag2 - k2*u_x_m_mag2*u_x_m_mag2;
        }
    }

}