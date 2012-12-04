cimport numpy as np

cdef extern from "clib.h":
	void compute_uniform_exch(double *spin,double *field,double J,double dx,double dy,double dz,int nx,int ny, int nz)
	void compute_anis(double *spin,double *field,double Kx,double Ky,double Kz,int nx,int ny, int nz)


def compute_uniform_exchange(
	np.ndarray[double, ndim=2, mode="c"] spin,
		np.ndarray[double, ndim=2, mode="c"] field,
		J,
		dx,dy,dz,
		nx,ny,nz
		):
	compute_uniform_exch(&spin[0,0],&field[0,0],J,dx,dy,dz,nx,ny,nz)    

def compute_anisotropy(
		np.ndarray[double, ndim=2, mode="c"] spin,
		np.ndarray[double, ndim=2, mode="c"] field,
		Kx,Ky,Kz,
		nx,ny,nz
		):
	compute_anis(&spin[0,0],&field[0,0],Kx,Ky,Kz,nx,ny,nz)    
