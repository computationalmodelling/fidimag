import numpy as np
cimport numpy as np

cdef extern from "clib.h":
	void compute_uniform_exch(double * spin, double * field, double J, double dx, double dy, double dz, int nx, int ny, int nz)
	void compute_anis(double * spin, double * field, double Dx, double Dy, double Dz, int nxyz)
	void llg_rhs(double * dm_dt, double * spin, double * h, double *alpha, double gamma, int nxyz)
	
	
	
	# used for demag
	ctypedef struct fft_demag_plan:
		pass
	
	fft_demag_plan * create_plan()
	void finalize_plan(fft_demag_plan * plan)
	void init_plan(fft_demag_plan * plan, double mu_s, double dx, double dy, double dz, int nx,
		int ny, int nz)
	void compute_fields(fft_demag_plan * plan, double * spin, double * field)
	void exact_compute(fft_demag_plan * plan, double * spin, double * field)
	
	
	# used for sllg
	ctypedef struct ode_solver:
		pass
	
	ode_solver *create_ode_plan()
	void init_solver(ode_solver *s, double mu_s, int nxyz, double dt, double gamma)
	void finalize_ode_plan(ode_solver *plan)
	void run_step1(ode_solver *s, double *m, double *h, double *m_pred, double *T, double *alpha)
	void run_step2(ode_solver *s, double *m_pred, double *h, double *m, double *T, double *alpha)
	
	

def compute_uniform_exchange(np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            J,
                            dx, dy, dz,
                            nx, ny, nz):
	
	compute_uniform_exch(& spin[0], & field[0], J, dx, dy, dz, nx, ny, nz)    

def compute_anisotropy(np.ndarray[double, ndim=1, mode="c"] spin,
                        np.ndarray[double, ndim=1, mode="c"] field,
                        Kx, Ky, Kz, nxyz):
	compute_anis(& spin[0], & field[0], Kx, Ky, Kz, nxyz)    


def compute_llg_rhs(np.ndarray[double, ndim=1, mode="c"] dm_dt,
    		np.ndarray[double, ndim=1, mode="c"] spin,
                np.ndarray[double, ndim=1, mode="c"] field,
                np.ndarray[double, ndim=1, mode="c"] alpha,
                gamma, nxyz):
	llg_rhs(&dm_dt[0], &spin[0], &field[0], &alpha[0], gamma, nxyz)
	
	

cdef class FFTDemag(object):
	cdef fft_demag_plan * _c_plan

	def __cinit__(self, mu_s, dx, dy, dz, nx, ny, nz):
		self._c_plan = create_plan()
		if self._c_plan is NULL:
			raise MemoryError()
		init_plan(self._c_plan, mu_s, dx, dy, dz, nx, ny, nz)
		
	def free(self):
		self.__dealloc__()
	
	def __dealloc__(self):
		if self._c_plan is not NULL:
			finalize_plan(self._c_plan)
			self._c_plan = NULL
			
			
	def compute_field(self,
					np.ndarray[double, ndim=1, mode="c"] spin,
					np.ndarray[double, ndim=1, mode="c"] field):
		compute_fields(self._c_plan, & spin[0], & field[0])
		
	def compute_exact(self,
					np.ndarray[double, ndim=1, mode="c"] spin,
                	np.ndarray[double, ndim=1, mode="c"] field):
		exact_compute(self._c_plan, & spin[0], & field[0])


		
cdef class RK2S(object):
	cdef ode_solver * _c_plan
	cdef double dt
	cdef update_fun
	cdef np.ndarray spin
	cdef np.ndarray pred_m
	cdef np.ndarray field
	cdef np.ndarray T
	cdef np.ndarray alpha
	
	cdef public double t
	cdef public np.ndarray y
	
	def __cinit__(self,mu_s,dt,nxyz,gamma,
				np.ndarray[double, ndim=1, mode="c"] alpha,
				np.ndarray[double, ndim=1, mode="c"] spin,
				np.ndarray[double, ndim=1, mode="c"] field,
				np.ndarray[double, ndim=1, mode="c"] T,
				update_fun):
			
		self.t = 0
		self.dt = dt
		self.y = spin
		self.update_fun = update_fun
		self.field = field
		self.T = T
		self.alpha= alpha
		self.pred_m=np.zeros(3*nxyz,dtype=np.float)
				
		self._c_plan = create_ode_plan()
		if self._c_plan is NULL:
			raise MemoryError()
		
		init_solver(self._c_plan,mu_s,nxyz,dt,gamma)
		
	def __dealloc__(self):
		if self._c_plan is not NULL:
			finalize_ode_plan(self._c_plan)
			self._c_plan = NULL
		
	
	def set_initial_value(self,np.ndarray[double, ndim=1, mode="c"] spin, t):
		self.t = t
		self.spin = spin
		

	def successful(self):
		#print self.spin
		return True
	
	def run_step(self, t):
		cdef np.ndarray[double, ndim=1, mode="c"] spin=self.spin
		cdef np.ndarray[double, ndim=1, mode="c"] field=self.field
		cdef np.ndarray[double, ndim=1, mode="c"] pred_m=self.pred_m
		cdef np.ndarray[double, ndim=1, mode="c"] T=self.T
		cdef np.ndarray[double, ndim=1, mode="c"] alpha=self.alpha
		
		#print "from cython1", self.spin,self.field,self.pred_m
		self.update_fun(self.spin)
		run_step1(self._c_plan,&spin[0],&field[0],&pred_m[0],&T[0],&alpha[0])
		#print "from cython2", self.spin,self.field,self.pred_m
		
		self.update_fun(self.pred_m)
		run_step2(self._c_plan,&pred_m[0],&field[0],&spin[0],&T[0],&alpha[0])
		self.t=t
					
	def integrate(self, t):
		while self.t<t:
			self.run_step(self.t+self.dt)
	
	
		
