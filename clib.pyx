cimport numpy as np

cdef extern from "clib.h":
	void compute_uniform_exch(double * spin, double * field, double J, double dx, double dy, double dz, int nx, int ny, int nz)
	void compute_anis(double *spin, double *field, double Dx, double Dy, double Dz,int nxyz)
	void llg_rhs(double * dm_dt, double * spin, double * h, double gamma, double alpha, double mu_s, int nxyz, double c)
	void compute_tensors(double *Nxx, double *Nyy, double *Nzz, double *Nxy,
		double *Nxz, double *Nyz, double dx, double dy, double dz, int nx,
		int ny, int nz)

def compute_uniform_exchange(
							np.ndarray[double, ndim=1, mode="c"] spin,
							np.ndarray[double, ndim=1, mode="c"] field,
							J,
                			dx, dy, dz,
                   			nx, ny, nz):
	compute_uniform_exch(& spin[0], & field[0], J, dx, dy, dz, nx, ny, nz)    

def compute_anisotropy(
					np.ndarray[double, ndim=1, mode="c"] spin,
                	np.ndarray[double, ndim=1, mode="c"] field,
                 	Kx, Ky, Kz,nxyz):
	compute_anis(& spin[0], & field[0], Kx, Ky, Kz, nxyz)    


def compute_llg_rhs(
				np.ndarray[double, ndim=1, mode="c"] dm_dt,
				np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                gamma, alpha, mu_s, nxyz, c):
	llg_rhs(& dm_dt[0], & spin[0], & field[0], gamma, alpha, mu_s, nxyz, c)
	
	
def compute_all_tensors(
				np.ndarray[double, ndim=1, mode="c"] Nxx,
				np.ndarray[double, ndim=1, mode="c"] Nyy,
				np.ndarray[double, ndim=1, mode="c"] Nzz,
				np.ndarray[double, ndim=1, mode="c"] Nxy,
				np.ndarray[double, ndim=1, mode="c"] Nxz,
				np.ndarray[double, ndim=1, mode="c"] Nyz,
                dx,dy,dz,nx,ny,nz):
	compute_tensors(&Nxx[0], &Nyy[0], &Nzz[0],&Nxy[0], &Nxz[0], &Nyz[0],dx,dy,dz,nx,ny,nz)
