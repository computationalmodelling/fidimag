
#include "clib.h"

void compute_anis(double *spin,double *field,double Kx,double Ky,double Kz,int nx,int ny, int nz){
  int nxy=nx*ny;
  int nxyz=nxy*nz;
  int i,j,k;
  int id;

  for (i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      for(k=0;k<nz;k++){
	id=nxy*k+ny*j+i;	
	field[id]=Kx*spin[id];
	id+=nxyz;
	field[id]=Ky*spin[id];
	id+=nxyz;
	field[id]=Kz*spin[id];
      }
    }
  }
  
}
