cimport numpy as np

cdef extern from "clib.h":
	void compute_uniform_exch(double *spin,double *field,double J,double dx,double dy,double dz,int nx,int ny, int nz)
	void compute_anis(double *spin,double *field,double Kx,double Ky,double Kz,int nx,int ny, int nz)
	void llg_rhs(double *dm_dt,double *spin,double *h,double gamma,double alpha, double mu_s, int nxyz)

def compute_uniform_exchange(
							np.ndarray[double, ndim=1, mode="c"] spin,
							np.ndarray[double, ndim=1, mode="c"] field,
							J,
                			dx,dy,dz,
                   			nx,ny,nz):
	compute_uniform_exch(&spin[0],&field[0],J,dx,dy,dz,nx,ny,nz)    

def compute_anisotropy(
					np.ndarray[double, ndim=1, mode="c"] spin,
                	np.ndarray[double, ndim=1, mode="c"] field,
                 	Kx,Ky,Kz,
                  	nx,ny,nz):
	compute_anis(&spin[0],&field[0],Kx,Ky,Kz,nx,ny,nz)    


def compute_llg_rhs(
				np.ndarray[double, ndim=1, mode="c"] dm_dt,
				np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                gamma,alpha,mu_s,nxyz):
	llg_rhs(&dm_dt[0],&spin[0],&field[0],gamma,alpha,mu_s,nxyz)