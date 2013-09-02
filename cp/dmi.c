#include "clib.h"

inline void cross_incr(double p[3], double q[3], double res[3]){
    res[0]+= -p[2]*q[1] + p[1]*q[2];
    res[1]+= p[2]*q[0] - p[0]*q[2];
    res[2]+= -p[1]*q[0] + p[0]*q[1];
}

void compute_dmi(double *spin, double *field, double Dx, double Dy, double Dz,
                          int nx, int ny, int nz) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int index, id;
	double f[3];
    double S[3];
	double D[3]={Dx,Dy,Dz};
    
    for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				index = nyz * i + nz * j + k;
				f[0] = 0;
				f[1] = 0;
				f[2] = 0;

				if (k > 0) {
					id = index - 1;
                    S[0] = spin[id];
                    S[1] = spin[id + n1];
                    S[2] = spin[id + n2];
					cross_incr(D,S,f);
				}

				if (j > 0) {
					id = index - nz;
                    S[0] = spin[id];
                    S[1] = spin[id + n1];
                    S[2] = spin[id + n2];
					cross_incr(D,S,f);
				}

				if (i > 0) {
					id = index - nyz;
                    S[0] = spin[id];
                    S[1] = spin[id + n1];
                    S[2] = spin[id + n2];
					cross_incr(D,S,f);
				}

				if (i < nx - 1) {
					id = index + nyz;
                    S[0] = spin[id];
                    S[1] = spin[id + n1];
                    S[2] = spin[id + n2];
					cross_incr(D,S,f);
				}

				if (j < ny - 1) {
					id = index + nz;
                    S[0] = spin[id];
                    S[1] = spin[id + n1];
                    S[2] = spin[id + n2];
					cross_incr(D,S,f);
				}

				if (k < nz - 1) {
					id = index + 1;
                    S[0] = spin[id];
                    S[1] = spin[id + n1];
                    S[2] = spin[id + n2];
					cross_incr(D,S,f);
				}

				field[index] = f[0];
				field[index + n1] = f[1];
				field[index + n2] = f[2];

			}
		}
	}

}

inline double single_energy(double Dv[3], double Si[3], double Sj[3]){
    double tx = Dv[2]*(-Si[1]*Sj[0]+Si[0]*Sj[1]);
    double ty = Dv[1]*(Si[2]*Sj[0]-Si[0]*Sj[2]);
    double tz = Dv[0]*(-Si[2]*Sj[1]+Si[1]*Sj[2]);
    return tx+ty+tz;
}

double compute_dmi_eny(double *spin, double Dx, double Dy, double Dz,
                          int nx, int ny, int nz) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int index, id;
	
	double energy = 0;
    double Dv[3]={Dx,Dy,Dz};
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
					energy += single_energy(Dv,S_i,S_j);
				}

				if (j < ny - 1) {
					id = index + nz;
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy(Dv,S_i,S_j);
				}

				if (k < nz - 1) {
					id = index + 1;
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy(Dv,S_i,S_j);
				}

			}
		}
	}

	energy = -energy;

	return energy;

}
