import numpy as np
cimport numpy as np

from cvode cimport *
from libc.string cimport  memcmp, memcpy, memmove

cdef extern from "cvode/cvode_spgmr.h":
    int CVSpgmr(void *cvode_mem, int pretype, int max1)

cdef extern from "cvode/cvode_diag.h":
    int CVDiag(void *cvode_mem)

cdef extern from "cvode/cvode_spils.h":
    int CVSpilsSetPrecType(void *cvode_mem, int pretype)
    int CVSpilsSetGSType(void *cvode_mem, int gstype)

cdef extern from "clib.h":
    void compute_uniform_exch(double * spin, double * field, double J, double dx, double dy, double dz, int nx, int ny, int nz)
    void compute_anis(double * spin, double * field, double Dx, double Dy, double Dz, int nxyz)
    void llg_rhs(double * dm_dt, double * spin, double * h, double *alpha, double gamma, int nxyz)
    int check_array(double *a, double *b, int n)


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


    def compute_field(self,np.ndarray[double, ndim=1, mode="c"] spin,
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

        self.update_fun = update_fun
        self.field = field
        self.T = T
        self.alpha= alpha
        self.pred_m=np.zeros(3*nxyz,dtype=np.float)
        self.y=np.zeros(3*nxyz,dtype=np.float)


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
        self.y[:] = spin[:]

    def successful(self):
        #print self.spin
        return True

    def run_step(self, t):
        cdef np.ndarray[double, ndim=1, mode="c"] y=self.y
        cdef np.ndarray[double, ndim=1, mode="c"] field=self.field
        cdef np.ndarray[double, ndim=1, mode="c"] pred_m=self.pred_m
        cdef np.ndarray[double, ndim=1, mode="c"] T=self.T
        cdef np.ndarray[double, ndim=1, mode="c"] alpha=self.alpha

        #print "from cython1", self.spin,self.field,self.pred_m
        self.update_fun(self.y)
        run_step1(self._c_plan,&y[0],&field[0],&pred_m[0],&T[0],&alpha[0])

        self.update_fun(self.pred_m)
        run_step2(self._c_plan,&pred_m[0],&field[0],&y[0],&T[0],&alpha[0])
        self.t = t

    def integrate(self, t):
        while self.t<t:
            self.run_step(self.t+self.dt)



cdef struct cv_userdata:
    int nxyz
    double mu_s
    double gamma
    void *alpha
    void *heff
    void *rhs_fun

cdef inline np.ndarray nv2arr(N_Vector v):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef double* v_data = (<N_VectorContent_Serial>v.content).data
    cdef np.ndarray[double, ndim=1, mode='c'] x=np.empty(n)
    memcpy(x.data, v_data, n*sizeof(double))
    return x


cdef int cv_rhs(realtype t, N_Vector yv, N_Vector yvdot, void* user_data):

    cdef cv_userdata *ud = <cv_userdata *>user_data

    cdef double *dm_dt = (<N_VectorContent_Serial>yvdot.content).data
    cdef int i, size = (<N_VectorContent_Serial>yvdot.content).length

    cdef np.ndarray[double, ndim=1, mode="c"] y=nv2arr(yv)
    cdef np.ndarray[double, ndim=1, mode="c"] heff=<np.ndarray>ud.heff
    cdef np.ndarray[double, ndim=1, mode="c"] alpha=<np.ndarray>ud.alpha
    
    #weird thing is cython gets error without this
    try:
        (<object>ud.rhs_fun)(y)
    except:
        pass

    llg_rhs(dm_dt, &y[0], &heff[0], &alpha[0], ud.gamma, ud.nxyz)

    return 0

cdef class CvodeLLG(object):
    
    cdef np.ndarray field
    cdef np.ndarray alpha

    cdef public double t
    cdef public np.ndarray y
    cdef N_Vector u_y

    cdef void *cvode_mem
    cdef void *rhs_fun
    cdef callback_fun
    cdef cv_userdata user_data
    cdef int MODIFIED_GS

    def __cinit__(self,mu_s,nxyz,gamma,
                            np.ndarray[double, ndim=1, mode="c"] alpha,
                            np.ndarray[double, ndim=1, mode="c"] spin,
                            np.ndarray[double, ndim=1, mode="c"] field,
                            callback_fun):

        self.t = 0

        self.callback_fun = callback_fun
        self.field = field
        self.alpha = alpha

        self.y = np.zeros(3*nxyz,dtype=np.float)
        cdef np.ndarray[double, ndim=1, mode="c"] y=self.y

        self.u_y = N_VMake_Serial(3*nxyz,&y[0])

        self.rhs_fun = <void *>cv_rhs

        self.user_data = cv_userdata(nxyz,mu_s,gamma,
                                     <void*>self.alpha,
                                     <void*>self.field,
                                     <void*>self.y,
                                     <void*>self.callback_fun)

        self.MODIFIED_GS = 1
        
        self.init_ode()

    def init_ode(self):
        self.cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);


    def set_initial_value(self,np.ndarray[double, ndim=1, mode="c"] spin, t):
        self.t = t
        self.y[:] = spin[:]

        flag = CVodeSetUserData(self.cvode_mem, <void*>&self.user_data);
        self.check_flag(flag,"CVodeSetUserData")

        flag = CVodeInit(self.cvode_mem, <CVRhsFn>self.rhs_fun, t, self.u_y)
        self.check_flag(flag,"CVodeInit")
        #print "m0:",nv2arr(self.u_y)

        flag = CVodeSStolerances(self.cvode_mem, 1e-6, 1e-6)

        #flag = CVSpgmr(self.cvode_mem, 1, 300);
        #flag = CVSpilsSetGSType(self.cvode_mem, 1);
        flag = CVDiag(self.cvode_mem)

    cpdef int run_until(self, double tf):
        cdef int flag
        cdef double tret
        flag = CVodeStep(self.cvode_mem, tf, self.u_y, &tret, CV_NORMAL)

        return 0


    def check_flag(self, flag, fun_name):
        if flag<0:
            raise Exception("Run %s failed!"%fun_name)
