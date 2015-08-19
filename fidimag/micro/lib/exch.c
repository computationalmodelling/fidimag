#include "micro_clib.h"

void compute_exch_field_micro(double *m, double *field, double *energy,
		double *Ms_inv, double A, double dx, double dy, double dz, int nx,
		int ny, int nz, int xperiodic, int yperiodic, int zperiodic) {

    int nxy = nx*ny;
    int nxyz = nxy*nz;

	double ax = 2 * A / (dx * dx), ay = 2 * A / (dy * dy), az = 2 * A / (dz * dz);

	#pragma omp parallel for
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
	       for (int i = 0; i < nx; i++) { 

				int index = nxy * k + nx * j + i;
				int indexm = 3*index;

                if (Ms_inv[index] == 0.0){
                    field[3*index] = 0;
                    field[3*index+1] = 0;
                    field[3*index+2] = 0;
                    energy[index] = 0;
                    continue;
                }

                int id=0, idm=0;
				double fx = 0, fy = 0, fz = 0;

                if (k > 0 || zperiodic ) {
                    id = index - nxy;
                    if (k == 0){
                        id += nxyz;
                    }
                    idm = 3*id;
					//in general this is unnecessary for normal LLG equation,
					//but for LLbar this is very important.
					if (Ms_inv[id] > 0) {
						fx += az * (m[idm] - m[indexm]);
						fy += az * (m[idm+1] - m[indexm+1]);
						fz += az * (m[idm+2] - m[indexm+2]);
					}
				}

                if (j > 0 || yperiodic) {
                    id = index - nx;
                    if (j==0) {
                        id += nxy;
                    }
                    idm = 3*id;
					if (Ms_inv[id] > 0) {
						fx += ay * (m[idm] - m[indexm]);
						fy += ay * (m[idm+1] - m[indexm+1]);
						fz += ay * (m[idm+2] - m[indexm+2]);
					}
				}

                if (i > 0 || xperiodic) {
                    id = index - 1;
                    if (i==0) {
                        id += nx;
                    }
                    idm = 3*id;
					if (Ms_inv[id] > 0) {
						fx += ax * (m[idm] - m[indexm]);
						fy += ax * (m[idm+1] - m[indexm+1]);
						fz += ax * (m[idm+2] - m[indexm+2]);
					}
				}

                if (i < nx - 1 || xperiodic) {
                    id = index + 1;
                    if (i == nx-1){
                        id -= nx;
                    }
                    idm = 3*id;
					if (Ms_inv[id] > 0) {
						fx += ax * (m[idm] - m[indexm]);
						fy += ax * (m[idm+1] - m[indexm+1]);
						fz += ax * (m[idm+2] - m[indexm+2]);
					}
				}

                if (j < ny - 1 || yperiodic) {
                    id = index + nx;
                    if (j == ny-1){
                        id -= nxy;
                    }
                    idm = 3*id;
					if (Ms_inv[id] > 0) {
						fx += ay * (m[idm] - m[indexm]);
						fy += ay * (m[idm+1] - m[indexm+1]);
						fz += ay * (m[idm+2] - m[indexm+2]);
					}
				}

                if (k < nz - 1 || zperiodic ) {
                    id = index + nxy;
                    if (k == nz-1){
                        k -= nxyz;
                    }
                    idm = 3*id;
					if (Ms_inv[id] > 0) {
						fx += az * (m[idm] - m[indexm]);
						fy += az * (m[idm+1] - m[indexm+1]);
						fz += az * (m[idm+2] - m[indexm+2]);
					}
				}

                energy[index] = -0.5*(fx*m[indexm]+fy*m[indexm+1]+fz*m[indexm+2]);
                
                // Note: both here and Dx don't have factor of 2.
                field[indexm] = fx*Ms_inv[index]*MU0_INV;
                field[indexm + 1] = fy*Ms_inv[index]*MU0_INV;
                field[indexm + 2] = fz*Ms_inv[index]*MU0_INV;

			}
		}
	}
}

