#include "clib.h"

void compute_uniform_exch(double *spin, double *field, double J, double dx,
                double dy, double dz, int nx, int ny, int nz) {

        int nxy = nx * ny;
        int nxyz = nxy * nz;
        int i, j, k;
        int index, id;
        double tmp[3];
        for (k = 0; k < nz; k++){
                for (j = 0; j < ny; j++) {
                	for (i = 0; i < nx; i++) {
                                index = nxy * k + nx * j + i;
                                tmp[0] = 0;
                                tmp[1] = 0;
                                tmp[2] = 0;

                                if (k > 0) {
                                        id = index - nxy;
                                        tmp[0] += J * spin[id];
                                        id += nxyz;
                                        tmp[1] += J * spin[id];
                                        id += nxyz;
                                        tmp[2] += J * spin[id];
                                }

                                if (j > 0) {
                                        id = index - nx;
                                        tmp[0] += J * spin[id];
                                        id += nxyz;
                                        tmp[1] += J * spin[id];
                                        id += nxyz;
                                        tmp[2] += J * spin[id];
                                }

                                if (i > 0) {
                                        id = index - 1;
                                        tmp[0] += J * spin[id];
                                        id += nxyz;
                                        tmp[1] += J * spin[id];
                                        id += nxyz;
                                        tmp[2] += J * spin[id];
                                }

                                if (i < nx - 1) {
                                        id = index + 1;
                                        tmp[0] += J * spin[id];
                                        id += nxyz;
                                        tmp[1] += J * spin[id];
                                        id += nxyz;
                                        tmp[2] += J * spin[id];
                                }

                                if (j < ny - 1) {
                                        id = index + nx;
                                        tmp[0] += J * spin[id];
                                        id += nxyz;
                                        tmp[1] += J * spin[id];
                                        id += nxyz;
                                        tmp[2] += J * spin[id];
                                }

                                if (k < nz - 1) {
                                        id = index + nxy;
                                        tmp[0] += J * spin[id];
                                        id += nxyz;
                                        tmp[1] += J * spin[id];
                                        id += nxyz;
                                        tmp[2] += J * spin[id];
                                }

                                field[index] = tmp[0];
                                index += nxyz;
                                field[index] = tmp[1];
                                index += nxyz;
                                field[index] = tmp[2];

                        }
                }
        }


}





double compute_exch_energy(double *spin, double J, int nx, int ny, int nz) {

        int nxy = nx * ny;
        int nxyz1 = nxy * nz, nxyz2 = 2*nxyz1;
        int i, j, k;
        int index, id;
        double Sx,Sy,Sz;
        double energy=0;

        for (k = 0; k < nz; k++){
                for (j = 0; j < ny; j++) {
                	for (i = 0; i < nx; i++) {
                		index = nxy * k + nx * j + i;
                		Sx = spin[index];
                		Sy = spin[index+nxyz1];
                		Sz = spin[index+nxyz2];

                		if (k > 0) {
                			id = index - nxy;
                			energy += J * Sx*spin[id];
                			energy += J * Sy*spin[id+nxyz1];
                			energy += J * Sz*spin[id+nxyz2];
                		}

                		if (j > 0) {
                			id = index - nx;
                			energy += J * Sx*spin[id];
                			energy += J * Sy*spin[id+nxyz1];
                			energy += J * Sz*spin[id+nxyz2];
                		}

                		if (i > 0) {
                			id = index - 1;
                			energy += J * Sx*spin[id];
                			energy += J * Sy*spin[id+nxyz1];
                			energy += J * Sz*spin[id+nxyz2];
                		}

                		if (i < nx - 1) {
                			id = index + 1;
                			energy += J * Sx*spin[id];
                			energy += J * Sy*spin[id+nxyz1];
                			energy += J * Sz*spin[id+nxyz2];
                		}

                		if (j < ny - 1) {
                			id = index + nx;
                			energy += J * Sx*spin[id];
                			energy += J * Sy*spin[id+nxyz1];
                			energy += J * Sz*spin[id+nxyz2];
                		}

                		if (k < nz - 1) {
                			id = index + nxy;
                			energy += J * Sx*spin[id];
                			energy += J * Sy*spin[id+nxyz1];
                			energy += J * Sz*spin[id+nxyz2];
                		}

                     }
                }
        }

        energy=-energy;
        //printf("energy=%g\n",energy);

        return energy;

}
