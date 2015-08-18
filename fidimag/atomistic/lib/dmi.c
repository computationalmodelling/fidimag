#include "clib.h"

void dmi_field_bulk(double *spin, double *field, double *energy, double D, double *ngbs, int nxyz) {
    
	#pragma omp parallel for
	for (int i = 0; i < nxyz; i++) {

		int id =0;
		int idv = 6*i; //index for the neighbours
		
		double fx = 0, fy=0, fz=0;

		if (ngbs[idv]>=0) { //x-1
			id = 3*ngbs[idv]; 
			fx += D*cross_x(-1,0,0,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(-1,0,0,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(-1,0,0,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+1]>=0) { //x+1
			id = 3*ngbs[idv+1]; 
			fx += D*cross_x(1,0,0,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(1,0,0,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(1,0,0,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+2]>=0) { //y-1
			id = 3*ngbs[idv+2]; 
			fx += D*cross_x(0,-1,0,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(0,-1,0,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(0,-1,0,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+3]>=0) { //y+1
			id = 3*ngbs[idv+3]; 
			fx += D*cross_x(0,1,0,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(0,1,0,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(0,1,0,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+4]>=0) { //z-1
			id = 3*ngbs[idv+4]; 
			fx += D*cross_x(0,0,-1,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(0,0,-1,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(0,0,-1,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+5]>=0) { //z+1
			id = 3*ngbs[idv+5]; 
			fx += D*cross_x(0,0,1,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(0,0,1,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(0,0,1,spin[id],spin[id+1],spin[id+2]);
		}

		field[3*i] = fx;
        field[3*i+1] = fy;
        field[3*i+2] = fz;
                
        energy[i] = -0.5*(fx*spin[3*i]+fy*spin[3*i+1]+fz*spin[3*i+2]);
    }

}

void dmi_field_interfacial_atomistic(double *m, double *field, double *energy,
    double D, int nx, int ny, int nz, int xperiodic, int yperiodic) {
    
  	int nyz = ny * nz;
  	int n1 = nx * nyz, n2 = 2 * n1;
    
  	double Dx = D, Dy = D;
    
  	#pragma omp parallel for
  	for (int i = 0; i < nx; i++) {
    	for (int j = 0; j < ny; j++) {
      		for (int k = 0; k < nz; k++) {
                
      		int index = nyz * i + nz * j + k;
                
      		double fx=0, fy=0, fz=0;
                
      		int id=0;
                
   	if (j > 0 || yperiodic) {
    	id = index - nz;
     	if (j==0) {id += nyz;}
     	//rij = (0,-1,0), rij \times z = (-1,0,0)
     	fx += Dy * cross_x(-1,0,0,m[id],m[id+n1],m[id+n2]);
     	fy += Dy * cross_y(-1,0,0,m[id],m[id+n1],m[id+n2]);
     	fz += Dy * cross_z(-1,0,0,m[id],m[id+n1],m[id+n2]);
  	}
                
  	if (i > 0 || xperiodic) {
    	id = index - nyz;
    	if (i==0) {id += n1;}
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

	field[index] = fx;
  	field[index + n1] = fy;
  	field[index + n2] = fz;
    
    //TODO: check whether the energy is correct or not. 
  	energy[index] = -0.5*(fx*m[index]+fy*m[index+n1]+fz*m[index+n2]);
                               
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
