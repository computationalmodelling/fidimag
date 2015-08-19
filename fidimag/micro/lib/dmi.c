#include "micro_clib.h"

void dmi_field_bulk(double *m, double *field, double *energy, double *Ms_inv,
                              double *D, double dx, double dy, double dz,
                              int nx, int ny, int nz, int xperiodic, int yperiodic, int zperiodic) {

    int nxy = nx*ny;
    int nxyz = nxy*nz;
    
	#pragma omp parallel for
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
	       for (int i = 0; i < nx; i++) { 
                
                int index = nxy * k + nx * j + i;

                double Dx = -D[index]/dx, Dy=-D[index]/dy, Dz=-D[index]/dz;
                
                double fx=0, fy=0, fz=0;
                
                int id = 0, idm=0; 
                
                if (Ms_inv[index] == 0.0){
                    field[3*index] = 0;
                    field[3*index+1] = 0;
                    field[3*index+2] = 0;
                    energy[index] = 0;
                    continue;
                }
                
                if (k > 0 || zperiodic ) {
                    id = index - nxy;
                    if (k == 0){
                        id += nxyz;
                    }
                    idm = 3*id;
                    fx += Dz * cross_x(0,0,-1,m[idm],m[idm+1],m[idm+2]);
                    fy += Dz * cross_y(0,0,-1,m[idm],m[idm+1],m[idm+2]);
                    fz += Dz * cross_z(0,0,-1,m[idm],m[idm+1],m[idm+2]);
                }
                
                if (j > 0 || yperiodic) {
                    id = index - nx;
                    if (j==0) {
                        id += nxy;
                    }
                    idm = 3*id;
                    fx += Dy * cross_x(0,-1,0,m[idm],m[idm+1],m[idm+2]);
                    fy += Dy * cross_y(0,-1,0,m[idm],m[idm+1],m[idm+2]);
                    fz += Dy * cross_z(0,-1,0,m[idm],m[idm+1],m[idm+2]);
                }
                
                if (i > 0 || xperiodic) {
                    id = index - 1;
                    if (i==0) {
                        id += nx;
                    }
                    idm = 3*id;
                    fx += Dx * cross_x(-1,0,0,m[idm],m[idm+1],m[idm+2]);
                    fy += Dx * cross_y(-1,0,0,m[idm],m[idm+1],m[idm+2]);
                    fz += Dx * cross_z(-1,0,0,m[idm],m[idm+1],m[idm+2]);
                }
                
                if (i < nx - 1 || xperiodic) {
                    id = index + 1;
                    if (i == nx-1){
                        id -= nx;
                    }
                    idm = 3*id;
                    fx += Dx * cross_x(1,0,0,m[idm],m[idm+1],m[idm+2]);
                    fy += Dx * cross_y(1,0,0,m[idm],m[idm+1],m[idm+2]);
                    fz += Dx * cross_z(1,0,0,m[idm],m[idm+1],m[idm+2]);
                }
                
                if (j < ny - 1 || yperiodic) {
                    id = index + nx;
                    if (j == ny-1){
                        id -= nxy;
                    }
                    idm = 3*id;
                    fx += Dy * cross_x(0,1,0,m[idm],m[idm+1],m[idm+2]);
                    fy += Dy * cross_y(0,1,0,m[idm],m[idm+1],m[idm+2]);
                    fz += Dy * cross_z(0,1,0,m[idm],m[idm+1],m[idm+2]);
                }
                
                if (k < nz - 1 || zperiodic ) {
                    id = index + nxy;
                    if (k == nz-1){
                        k -= nxyz;
                    }
                    idm = 3*id;
                    fx += Dz * cross_x(0,0,1,m[idm],m[idm+1],m[idm+2]);
                    fy += Dz * cross_y(0,0,1,m[idm],m[idm+1],m[idm+2]);
                    fz += Dz * cross_z(0,0,1,m[idm],m[idm+1],m[idm+2]);
                }
                
                energy[index] = -0.5*(fx*m[3*index]+fy*m[3*index+1]+fz*m[3*index+2]);
                
                // Note: both here and Dx don't have factor of 2.
                field[3*index] = fx*Ms_inv[index]*MU0_INV;
                field[3*index + 1] = fy*Ms_inv[index]*MU0_INV;
                field[3*index + 2] = fz*Ms_inv[index]*MU0_INV;
                
            }
        }
	}

}


void dmi_field_interfacial(double *m, double *field, double *energy, double *Ms_inv,
                    double *D, double dx, double dy, double dz,
                    int nx, int ny, int nz, int xperiodic, int yperiodic) {
    
    int nxy = nx*ny;
    
    #pragma omp parallel for
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
           for (int i = 0; i < nx; i++) { 
                
                int index = nxy * k + nx * j + i;
                
                double Dx = D[index]/dx, Dy=D[index]/dy;
                
                double fx=0, fy=0, fz=0;
                
                int id=0, idm=0;
                
                if (Ms_inv[index] == 0.0){
                    field[3*index] = 0;
                    field[3*index+1] = 0;
                    field[3*index+2] = 0;
                    energy[index] = 0;
                    continue;
                }


                if (j > 0 || yperiodic) {
                    id = index - nx;
                    if (j==0) {
                        id += nxy;
                    }
                    idm = 3*id;
                    //rij = (0,-1,0), rij \times z = (-1,0,0)
                    fx += Dy * cross_x(-1,0,0,m[idm],m[idm+1],m[idm+2]);
                    fy += Dy * cross_y(-1,0,0,m[idm],m[idm+1],m[idm+2]);
                    fz += Dy * cross_z(-1,0,0,m[idm],m[idm+1],m[idm+2]);
                }
                
                if (i > 0 || xperiodic) {
                    id = index - 1;
                    if (i==0) {
                        id += nx;
                    }
                    idm = 3*id;
                    //rij = (-1,0,0), rij \times z = (0,1,0)
                    fx += Dx * cross_x(0,1,0,m[idm],m[idm+1],m[idm+2]);
                    fy += Dx * cross_y(0,1,0,m[idm],m[idm+1],m[idm+2]);
                    fz += Dx * cross_z(0,1,0,m[idm],m[idm+1],m[idm+2]);
                }
                
                if (i < nx - 1 || xperiodic) {
                    id = index + 1;
                    if (i == nx-1){
                        id -= nx;
                    }
                    idm = 3*id;
                    //rij = (1,0,0), rij \times z = (0,-1,0)
                    fx += Dx * cross_x(0,-1,0,m[idm],m[idm+1],m[idm+2]);
                    fy += Dx * cross_y(0,-1,0,m[idm],m[idm+1],m[idm+2]);
                    fz += Dx * cross_z(0,-1,0,m[idm],m[idm+1],m[idm+2]);
                }
                
                if (j < ny - 1 || yperiodic) {
                    id = index + nx;
                    if (j == ny-1){
                        id -= nxy;
                    }
                    idm = 3*id;
                    //rij = (0,1,0), rij \times z = (1,0,0)
                    fx += Dy * cross_x(1,0,0,m[idm],m[idm+1],m[idm+2]);
                    fy += Dy * cross_y(1,0,0,m[idm],m[idm+1],m[idm+2]);
                    fz += Dy * cross_z(1,0,0,m[idm],m[idm+1],m[idm+2]);
                }

                
                energy[index] = -0.5*(fx*m[3*index]+fy*m[3*index+1]+fz*m[3*index+2]);
                
                // Note: both here and Dx don't have factor of 2.
                field[3*index] = fx*Ms_inv[index]*MU0_INV;
                field[3*index + 1] = fy*Ms_inv[index]*MU0_INV;
                field[3*index + 2] = fz*Ms_inv[index]*MU0_INV;
                
            }
        }
	}
    
}
