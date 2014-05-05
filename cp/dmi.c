#include "clib.h"

inline double cross_x(double a0, double a1, double a2, double b0, double b1, double b2) { return a1*b2 - a2*b1; }
inline double cross_y(double a0, double a1, double a2, double b0, double b1, double b2) { return a2*b0 - a0*b2; }
inline double cross_z(double a0, double a1, double a2, double b0, double b1, double b2) { return a0*b1 - a1*b0; }


void dmi_field(double *spin, double *field, double *energy, double D, int nx, int ny, int nz, int xperiodic, int yperiodic) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int index, id;
	double fx,fy,fz;
    
	#pragma omp parallel for private(i,j,k, index, id, fx, fy, fz)
	for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                
                index = nyz * i + nz * j + k;
                fx = 0;
                fy = 0;
                fz = 0;
                
                if (k > 0) {
                    id = index - 1;
                    fx += D * cross_x(0,0,-1,spin[id],spin[id+n1],spin[id+n2]);
                    fy += D * cross_y(0,0,-1,spin[id],spin[id+n1],spin[id+n2]);
                    fz += D * cross_z(0,0,-1,spin[id],spin[id+n1],spin[id+n2]);
                }
                
                if (j > 0 || yperiodic) {
                    id = index - nz;
                    if (j==0) {
                        id += nyz;
                    }
                    fx += D * cross_x(0,-1,0,spin[id],spin[id+n1],spin[id+n2]);
                    fy += D * cross_y(0,-1,0,spin[id],spin[id+n1],spin[id+n2]);
                    fz += D * cross_z(0,-1,0,spin[id],spin[id+n1],spin[id+n2]);
                }
                
                if (i > 0 || xperiodic) {
                    id = index - nyz;
                    if (i==0) {
                        id += n1;
                    }
                    fx += D * cross_x(-1,0,0,spin[id],spin[id+n1],spin[id+n2]);
                    fy += D * cross_y(-1,0,0,spin[id],spin[id+n1],spin[id+n2]);
                    fz += D * cross_z(-1,0,0,spin[id],spin[id+n1],spin[id+n2]);
                }
                
                if (i < nx - 1 || xperiodic) {
                    id = index + nyz;
                    if (i == nx-1){
                        id -= n1;
                    }
                    fx += D * cross_x(1,0,0,spin[id],spin[id+n1],spin[id+n2]);
                    fy += D * cross_y(1,0,0,spin[id],spin[id+n1],spin[id+n2]);
                    fz += D * cross_z(1,0,0,spin[id],spin[id+n1],spin[id+n2]);
                }
                
                if (j < ny - 1 || yperiodic) {
                    id = index + nz;
                    if (j == ny-1){
                        id -= nyz;
                    }
                    fx += D * cross_x(0,1,0,spin[id],spin[id+n1],spin[id+n2]);
                    fy += D * cross_y(0,1,0,spin[id],spin[id+n1],spin[id+n2]);
                    fz += D * cross_z(0,1,0,spin[id],spin[id+n1],spin[id+n2]);
                }
                
                if (k < nz - 1) {
                    id = index + 1;
                    fx += D * cross_x(0,0,1,spin[id],spin[id+n1],spin[id+n2]);
                    fy += D * cross_y(0,0,1,spin[id],spin[id+n1],spin[id+n2]);
                    fz += D * cross_z(0,0,1,spin[id],spin[id+n1],spin[id+n2]);
                }
                
                field[index] = fx;
                field[index + n1] = fy;
                field[index + n2] = fz;
                
                energy[index] = -(fx*spin[index]+fy*spin[index+n1]+fz*spin[index+n2]);

            }
        }
	}

}

inline double single_energy_x(double D, double Si[3], double Sj[3]){
    double tx = D*(-Si[2]*Sj[1]+Si[1]*Sj[2]);
    return tx;
}

inline double single_energy_y(double D, double Si[3], double Sj[3]){
    double ty = D*(Si[2]*Sj[0]-Si[0]*Sj[2]);
    return ty;
}

inline double single_energy_z(double D, double Si[3], double Sj[3]){
    double tz = D*(-Si[1]*Sj[0]+Si[0]*Sj[1]);
    return tz;
}

double dmi_energy(double *spin, double D, int nx, int ny, int nz, int xperiodic, int yperiodic) {

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

				if (i < nx - 1 || xperiodic) {
					id = index + nyz;
                    if (i == nx -1){
                        id -= n1;
                    }
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_x(D,S_i,S_j);
				}

				if (j < ny - 1 || yperiodic) {
					id = index + nz;
                    if (j == ny-1){
                        id -= nyz;
                    }
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_y(D,S_i,S_j);
				}

				if (k < nz - 1) {
					id = index + 1;
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_z(D,S_i,S_j);
				}

			}
		}
	}

	return energy;

}
