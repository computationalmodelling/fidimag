import numpy as np
cimport numpy as np  # import special compile-time information about numpy
np.import_array()  # don't remove or you'll segfault
from libc.string cimport memcpy


cdef extern from "../../atomistic/lib/clib.h":
    void normalise(double *m, int nxyz)


cdef struct cv_userdata:
    void *rhs_fun
    void *y
    void *dm_dt
    void *jvn_fun
    void *v  # is the vector by which the Jacobian must be multiplied
    void *jv  # the value of Jacobian times v


cdef inline copy_arr2nv(np.ndarray[realtype, ndim=1, mode='c'] np_x, N_Vector v):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef void* data_ptr=<void *>np_x.data
    memcpy((<N_VectorContent_Serial>v.content).data, data_ptr, n*sizeof(double))
    return 0


cdef inline copy_nv2arr(N_Vector v, np.ndarray[realtype, ndim=1, mode='c'] np_x):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef double* v_data = (<N_VectorContent_Serial>v.content).data
    memcpy(np_x.data, v_data, n*sizeof(realtype))
    return 0


cdef int cv_rhs(realtype t, N_Vector yv, N_Vector yvdot, void* user_data) except -1:
    cdef cv_userdata *ud = <cv_userdata *>user_data
    cdef np.ndarray[realtype, ndim=1, mode='c'] y_arr= <np.ndarray[realtype, ndim=1, mode='c']>ud.y
    cdef np.ndarray ydot_arr= <np.ndarray[realtype, ndim=1, mode='c']>ud.dm_dt

    copy_nv2arr(yv, y_arr)
    (<object>ud.rhs_fun)(t, y_arr, ydot_arr)
    copy_arr2nv(ydot_arr, yvdot)
    return 0


cdef int cv_jtimes(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp) except -1:
    cdef cv_userdata *ud = <cv_userdata *>user_data
    cdef np.ndarray[realtype, ndim=1, mode='c'] m_arr = <np.ndarray[realtype, ndim=1, mode='c']>ud.y
    cdef np.ndarray[realtype, ndim=1, mode='c'] mp_arr = <np.ndarray[realtype, ndim=1, mode='c']>ud.v
    cdef np.ndarray[realtype, ndim=1, mode='c'] ydot_arr = <np.ndarray[realtype, ndim=1, mode='c']>ud.dm_dt
    cdef np.ndarray[realtype, ndim=1, mode='c'] jv_arr = <np.ndarray[realtype, ndim=1, mode='c']>ud.jv

    copy_nv2arr(y, m_arr)
    copy_nv2arr(v, mp_arr)
    copy_nv2arr(fy, ydot_arr)
    (<object>ud.jvn_fun)(mp_arr, jv_arr, t, m_arr, ydot_arr)
    copy_arr2nv(jv_arr, Jv)
    return 0


cdef int psolve(realtype t, N_Vector y, N_Vector fy,
                N_Vector r, N_Vector z, realtype gamma, realtype delta, int lr,
                void *user_data, N_Vector tmp):
    z[:] = r
    return 0


cdef class CvodeSolver(object):
    cdef public double t
    cdef public np.ndarray y
    cdef double rtol, atol
    cdef np.ndarray spin
    cdef np.ndarray dm_dt
    cdef np.ndarray mp
    cdef np.ndarray Jmp
    cdef N_Vector u_y
    cdef void *cvode_mem
    cdef void *rhs_fun
    cdef void *jvn_fun
    cdef callback_fun
    cdef jtimes_fun
    cdef cv_userdata user_data
    cdef long int nsteps, nfevals, njevals
    cdef int max_num_steps
    cdef int has_jtimes

    def __cinit__(self, spin, callback_fun, jtimes_fun=None, rtol=1e-8, atol=1e-8):
        self.t = 0
        self.spin = spin
        self.dm_dt = np.copy(spin)
        self.y = np.copy(spin)
        self.mp = np.copy(spin)
        self.Jmp = np.copy(spin)

        self.callback_fun = callback_fun
        self.jtimes_fun=jtimes_fun
        self.jvn_fun = <void *>cv_jtimes
        self.rhs_fun = <void *>cv_rhs # wrapper for callback_fun (which is a Python function)

        self.has_jtimes = 0
        if jtimes_fun is not None:
            self.has_jtimes = 1

        self.user_data = cv_userdata(<void*>self.callback_fun,
                                     <void *>self.spin, <void *>self.dm_dt,
                                     <void *>self.jtimes_fun,
                                     <void *>self.mp,<void *>self.Jmp)

        self.cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
        flag = CVodeSetUserData(self.cvode_mem, <void*>&self.user_data);
        self.check_flag(flag,"CVodeSetUserData")
        self.set_initial_value(spin, self.t, 0)
        self.set_options(rtol, atol)

    def reset(self, np.ndarray[double, ndim=1, mode="c"] spin, t):
        copy_arr2nv(spin, self.u_y)
        CVodeReInit(self.cvode_mem, t, self.u_y)

    def set_initial_value(self, np.ndarray[double, ndim=1, mode="c"] spin, t, flag_m):
        self.t = t
        self.y[:] = spin[:]

        cdef np.ndarray[double, ndim=1, mode="c"] y = self.y
        self.u_y = N_VMake_Serial(y.size, &y[0])

        # We make sure that new memory for the magnetisation data is only allocated
        # once, by using the `flag_m` variable to signify allocation has already
        # taken place. If `flag_m` is 1, we use CvodeReInit instead of CvodeInit,
        # which reuses existing memory.
        if not flag_m:
            flag = CVodeInit(self.cvode_mem, <CVRhsFn>self.rhs_fun, t, self.u_y)
            self.check_flag(flag, "CVodeInit")
        else:
            flag = CVodeReInit(self.cvode_mem, t, self.u_y)

        if self.has_jtimes:  # we will use preconditioning with the Jacobian
            # CVSpgmr(cvode_mem, pretype, maxl) p. 27 of CVODE 2.7 manual
            flag = CVSpgmr(self.cvode_mem, PREC_LEFT, 300);
            # functions below in p. 37 CVODE 2.7 manual
            flag = CVSpilsSetJacTimesVecFn(self.cvode_mem, <CVSpilsJacTimesVecFn>self.jvn_fun)
            flag = CVSpilsSetPreconditioner(self.cvode_mem, <CVSpilsPrecSetupFn>self.psetup, <CVSpilsPrecSolveFn>psolve)
        else:
            flag = CVSpgmr(self.cvode_mem, PREC_NONE, 300);

    def set_options(self, rtol, atol, max_num_steps=100000, max_ord=None):
        self.rtol = rtol
        self.atol = atol
        self.max_num_steps = max_num_steps

        # Set tolerances
        flag = CVodeSStolerances(self.cvode_mem, self.rtol, self.atol)

        # maximum order of the linear multistep method
        # c.f. Suess et al. / Journal of Magnetism and Magnetic Materials 248 (2002) 298–311
        # where it speeds up time integration considerably (Fig. 3)
        if max_ord is not None:  # default is 5 for BDF method, 2 is a good alternative value to try
            flag = CVodeSetMaxOrd(self.cvode_mem, max_ord)

        # Set maximum number of iteration steps (?)
        flag = CVodeSetMaxNumSteps(self.cvode_mem, max_num_steps)
        CVodeReInit(self.cvode_mem, self.t, self.u_y)

    cpdef int run_until(self, double tf) except -1:
        cdef int flag
        cdef double tret

        flag = CVodeStep(self.cvode_mem, tf, self.u_y, &tret, CV_NORMAL)
        self.check_flag(flag,"CVodeStep")
        self.t = tret  # FIXME: or should this be tf?
        return 0

    cdef int psetup(self, realtype t, N_Vector y, N_Vector fy,
                     booleantype jok, booleantype *jcurPtr, realtype gamma,
                     void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
        if not jok:
            copy_nv2arr(y, self.y)
        return 0

    def check_flag(self, flag, fun_name):
        if flag < 0:
            raise Exception("Run %s failed!"%fun_name)

    def stat(self):
        CVodeGetNumSteps(self.cvode_mem, &self.nsteps);
        CVodeGetNumRhsEvals(self.cvode_mem, &self.nfevals);
        CVSpilsGetNumJtimesEvals(self.cvode_mem, &self.njevals)
        return self.nsteps, self.nfevals, self.njevals

    def get_current_step(self):
        cdef double step
        CVodeGetCurrentStep(self.cvode_mem, &step)
        return step

    def __repr__(self):
        s = []
        s.append("nsteps = %d," % self.nsteps)
        s.append("nfevals = %d," % self.nfevals)
        s.append("njevals = %d.\n" % self.njevals)

        return "(%s)" % ("\n".join(s))

    def __str__(self):
        return '%s%s' % (self.__class__.__name__, self.__repr__())

    def __dealloc__(self):
        self.user_data.rhs_fun = NULL
        self.user_data.y = NULL
        self.user_data.jvn_fun = NULL
        self.user_data.dm_dt = NULL
        self.user_data.v = NULL
        self.user_data.jv = NULL

        N_VDestroy_Serial(self.u_y)
        CVodeFree(&self.cvode_mem)
