#include "micro_clib.h"

void dmi_field_bulk(double *m, double *field, double *energy, double *Ms_inv,
                              double D, double dx, double dy, double dz,
                              int nx, int ny, int nz, int xperiodic, int yperiodic) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
    
    double Dx = -D/dx, Dy=-D/dy, Dz=-D/dz;
    
	#pragma omp parallel for
	for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                
                int index = nyz * i + nz * j + k;
                
                double fx=0, fy=0, fz=0;
                
                int id=0;
                
                if (Ms_inv[index] == 0.0){
                    field[index] = 0;
                    field[index + n1] = 0;
                    field[index + n2] = 0;
                    energy[index] = 0;
                    continue;
                }
                
                if (k > 0) {
                    id = index - 1;
                    fx += Dz * cross_x(0,0,-1,m[id],m[id+n1],m[id+n2]);
                    fy += Dz * cross_y(0,0,-1,m[id],m[id+n1],m[id+n2]);
                    fz += Dz * cross_z(0,0,-1,m[id],m[id+n1],m[id+n2]);
                }
                
                if (j > 0 || yperiodic) {
                    id = index - nz;
                    if (j==0) {
                        id += nyz;
                    }
                    fx += Dy * cross_x(0,-1,0,m[id],m[id+n1],m[id+n2]);
                    fy += Dy * cross_y(0,-1,0,m[id],m[id+n1],m[id+n2]);
                    fz += Dy * cross_z(0,-1,0,m[id],m[id+n1],m[id+n2]);
                }
                
                if (i > 0 || xperiodic) {
                    id = index - nyz;
                    if (i==0) {
                        id += n1;
                    }
                    fx += Dx * cross_x(-1,0,0,m[id],m[id+n1],m[id+n2]);
                    fy += Dx * cross_y(-1,0,0,m[id],m[id+n1],m[id+n2]);
                    fz += Dx * cross_z(-1,0,0,m[id],m[id+n1],m[id+n2]);
                }
                
                if (i < nx - 1 || xperiodic) {
                    id = index + nyz;
                    if (i == nx-1){
                        id -= n1;
                    }
                    fx += Dx * cross_x(1,0,0,m[id],m[id+n1],m[id+n2]);
                    fy += Dx * cross_y(1,0,0,m[id],m[id+n1],m[id+n2]);
                    fz += Dx * cross_z(1,0,0,m[id],m[id+n1],m[id+n2]);
                }
                
                if (j < ny - 1 || yperiodic) {
                    id = index + nz;
                    if (j == ny-1){
                        id -= nyz;
                    }
                    fx += Dy * cross_x(0,1,0,m[id],m[id+n1],m[id+n2]);
                    fy += Dy * cross_y(0,1,0,m[id],m[id+n1],m[id+n2]);
                    fz += Dy * cross_z(0,1,0,m[id],m[id+n1],m[id+n2]);
                }
                
                if (k < nz - 1) {
                    id = index + 1;
                    fx += Dz * cross_x(0,0,1,m[id],m[id+n1],m[id+n2]);
                    fy += Dz * cross_y(0,0,1,m[id],m[id+n1],m[id+n2]);
                    fz += Dz * cross_z(0,0,1,m[id],m[id+n1],m[id+n2]);
                }
                
                energy[index] = -0.5*(fx*m[index]+fy*m[index+n1]+fz*m[index+n2]);
                
                // Note: both here and Dx don't have factor of 2.
                field[index] = fx*Ms_inv[index]*MU0_INV;
                field[index + n1] = fy*Ms_inv[index]*MU0_INV;
                field[index + n2] = fz*Ms_inv[index]*MU0_INV;
                
            }
        }
	}

}


void dmi_field_interfacial(double *m, double *field, double *energy, double *Ms_inv,
                    double D, double dx, double dy, double dz,
                    int nx, int ny, int nz, int xperiodic, int yperiodic) {
    
	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
    
    double Dx = D/dx, Dy=D/dy;
    
    #pragma omp parallel for
	for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                
                int index = nyz * i + nz * j + k;
                
                double fx=0, fy=0, fz=0;
                
                int id=0;
                
                if (Ms_inv[index] == 0.0){
                    field[index] = 0;
                    field[index + n1] = 0;
                    field[index + n2] = 0;
                    energy[index] = 0;
                    continue;
                }
                
                if (j > 0 || yperiodic) {
                    id = index - nz;
                    if (j==0) {
                        id += nyz;
                    }
                    //rij = (0,-1,0), rij \times z = (-1,0,0)
                    fx += Dy * cross_x(-1,0,0,m[id],m[id+n1],m[id+n2]);
                    fy += Dy * cross_y(-1,0,0,m[id],m[id+n1],m[id+n2]);
                    fz += Dy * cross_z(-1,0,0,m[id],m[id+n1],m[id+n2]);
                }
                
                if (i > 0 || xperiodic) {
                    id = index - nyz;
                    if (i==0) {
                        id += n1;
                    }
                    //rij = (-1,0,0), rij \times z = (0,1,0)

                    fx += Dx * cross_x(0,1,0,m[id],m[id+n1],m[id+n2]);
                    fy += Dx * cross_y(0,1,0,m[id],m[id+n1],m[id+n2]);
                    fz += Dx * cross_z(0,1,0,m[id],m[id+n1],m[id+n2]);
                }
                
                if (i < nx - 1 || xperiodic) {
                    id = index + nyz;
                    if (i == nx-1){
                        id -= n1;
                    }
                    //rij = (1,0,0), rij \times z = (0,-1,0)
                    
                    fx += Dx * cross_x(0,-1,0,m[id],m[id+n1],m[id+n2]);
                    fy += Dx * cross_y(0,-1,0,m[id],m[id+n1],m[id+n2]);
                    fz += Dx * cross_z(0,-1,0,m[id],m[id+n1],m[id+n2]);
                }
                
                if (j < ny - 1 || yperiodic) {
                    id = index + nz;
                    if (j == ny-1){
                        id -= nyz;
                    }
                    //rij = (0,1,0), rij \times z = (1,0,0)
                    fx += Dy * cross_x(1,0,0,m[id],m[id+n1],m[id+n2]);
                    fy += Dy * cross_y(1,0,0,m[id],m[id+n1],m[id+n2]);
                    fz += Dy * cross_z(1,0,0,m[id],m[id+n1],m[id+n2]);
                }
                
                energy[index] = -0.5*(fx*m[index]+fy*m[index+n1]+fz*m[index+n2]);
                
                // Note: both here and Dx don't have factor of 2.
                field[index] = fx*Ms_inv[index]*MU0_INV;
                field[index + n1] = fy*Ms_inv[index]*MU0_INV;
                field[index + n2] = fz*Ms_inv[index]*MU0_INV;
                
            }
        }
	}
    
}
