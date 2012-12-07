cimport numpy as np

cdef extern from "clib.h":
	void compute_uniform_exch(double * spin, double * field, double J, double dx, double dy, double dz, int nx, int ny, int nz)
	void compute_anis(double *spin, double *field, double Dx, double Dy, double Dz,int nxyz)
	void llg_rhs(double * dm_dt, double * spin, double * h, double gamma, double alpha, double mu_s, int nxyz, double c)
	
	
	
	#used for demag
	ctypedef struct fft_demag_plan:
		pass
	
	fft_demag_plan *create_plan()
	void finalize_plan(fft_demag_plan *plan)
	void init_plan(fft_demag_plan *plan, double dx, double dy, double dz, int nx,
		int ny, int nz)
	void compute_fields(fft_demag_plan *plan, double *spin, double *field)
	void exact_compute(fft_demag_plan *plan, double *spin, double *field)


def compute_uniform_exchange(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            J,
                            dx, dy, dz,
                            nx, ny, nz):
	compute_uniform_exch(& spin[0], & field[0], J, dx, dy, dz, nx, ny, nz)    

def compute_anisotropy(np.ndarray[double, ndim=1, mode="c"] spin,
                        np.ndarray[double, ndim=1, mode="c"] field,
                 	Kx, Ky, Kz,nxyz):
	compute_anis(& spin[0], & field[0], Kx, Ky, Kz, nxyz)    


def compute_llg_rhs(
				np.ndarray[double, ndim=1, mode="c"] dm_dt,
				np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                gamma, alpha, mu_s, nxyz, c):
	llg_rhs(& dm_dt[0], & spin[0], & field[0], gamma, alpha, mu_s, nxyz, c)
	
	

cdef class FFTDemag:
	cdef fft_demag_plan * _c_plan

	def __cinit__(self,dx,dy,dz,nx,ny,nz):
		self._c_plan = create_plan()
		print 'from cython p='
		if self._c_plan is NULL:
			raise MemoryError()
		init_plan(self._c_plan,dx,dy,dz,nx,ny,nz)
		
	def free(self):
		self.__dealloc__()
	
	def __dealloc__(self):
		if self._c_plan is not NULL:
			finalize_plan(self._c_plan)
			self._c_plan = NULL
			
			
	def compute_field(self,
					np.ndarray[double, ndim=1, mode="c"] spin,
                	np.ndarray[double, ndim=1, mode="c"] field):
		compute_fields(self._c_plan,&spin[0],&field[0])
		
	def compute_exact(self,
					np.ndarray[double, ndim=1, mode="c"] spin,
                	np.ndarray[double, ndim=1, mode="c"] field):
		exact_compute(self._c_plan,&spin[0],&field[0])
		
