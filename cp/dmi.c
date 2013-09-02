#include "clib.h"

void compute_dmi(double *spin, double *field, double Dx, double Dy, double Dz,
                          int nx, int ny, int nz) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int index, id;
	double fx,fy,fz;
    double Sx,Sy,Sz;
    
    for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				index = nyz * i + nz * j + k;
				fx = 0;
				fy = 0;
				fz = 0;

				if (k > 0) {
					id = index - 1;
                    Sx = spin[id];
                    Sy = spin[id + n1];
                    Sz = spin[id + n2];
                    fx += -Dz*Sy;
                    fy += Dz*Sx;
				}

				if (j > 0) {
					id = index - nz;
                    Sx = spin[id];
                    Sy = spin[id + n1];
                    Sz = spin[id + n2];
                    fx += Dy*Sz;
                    fz += -Dy*Sx;
				}

				if (i > 0) {
					id = index - nyz;
                    Sx = spin[id];
                    Sy = spin[id + n1];
                    Sz = spin[id + n2];
                    fy += -Dx*Sz;
                    fz += Dx*Sy;
				}

				if (i < nx - 1) {
					id = index + nyz;
                    Sx = spin[id];
                    Sy = spin[id + n1];
                    Sz = spin[id + n2];
                    fy += Dx*Sz;
                    fz += -Dx*Sy;
				}

				if (j < ny - 1) {
					id = index + nz;
                    Sx = spin[id];
                    Sy = spin[id + n1];
                    Sz = spin[id + n2];
                    fx += -Dy*Sz;
                    fz += Dy*Sx;
				}

				if (k < nz - 1) {
					id = index + 1;
                    Sx = spin[id];
                    Sy = spin[id + n1];
                    Sz = spin[id + n2];
                    fx += Dz*Sy;
                    fy += -Dz*Sx;
				}

				field[index] = fx;
				field[index + n1] = fy;
				field[index + n2] = fz;

			}
		}
	}

}

inline double single_energy_x(double Dx, double Si[3], double Sj[3]){
    double tx = Dx*(-Si[2]*Sj[1]+Si[1]*Sj[2]);
    return tx;
}

inline double single_energy_y(double Dy, double Si[3], double Sj[3]){
    double ty = Dy*(Si[2]*Sj[0]-Si[0]*Sj[2]);
    return ty;
}

inline double single_energy_z(double Dz, double Si[3], double Sj[3]){
    double tz = Dz*(-Si[1]*Sj[0]+Si[0]*Sj[1]);
    return tz;
}

double compute_dmi_eny(double *spin, double Dx, double Dy, double Dz,
                          int nx, int ny, int nz) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int index, id;
	
	double energy = 0;
    
    double S_i[3],S_j[3];

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				index = nyz * i + nz * j + k;
				S_i[0] = spin[index];
				S_i[1] = spin[index + n1];
				S_i[2] = spin[index + n2];

				if (i < nx - 1) {
					id = index + nyz;
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_x(Dx,S_i,S_j);
				}

				if (j < ny - 1) {
					id = index + nz;
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_y(Dy,S_i,S_j);
				}

				if (k < nz - 1) {
					id = index + 1;
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_z(Dz,S_i,S_j);
				}

			}
		}
	}

	return energy;

}
