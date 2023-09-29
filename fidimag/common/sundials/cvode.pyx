import numpy as np
import os
cimport numpy as np  # import special compile-time information about numpy
cimport openmp
np.import_array()  # don't remove or you'll segfault
from libc.string cimport memcpy
import sys

cdef extern from "sundials/sundials_context.h":
    cdef struct _SUNContext
    ctypedef _SUNContext *SUNContext
    int SUNContext_Create(void *comm, SUNContext *ctx)
    int SUNContext_Free(SUNContext *ctx)

cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype
    ctypedef bint booleantype

cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void *content

    ctypedef _generic_N_Vector *N_Vector
    N_Vector N_VNew_Serial(long int vec_length, SUNContext sunctx)
    N_Vector N_VNew_OpenMP(long int vec_length, int num_threads, SUNContext sunctx)
    void N_VDestroy_Serial(N_Vector v)
    void N_VDestroy_OpenMP(N_Vector v)
    void N_VPrint_Serial(N_Vector v)
    void N_VPrint_OpenMP(N_Vector v)

cdef extern from "nvector/nvector_serial.h":
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data, SUNContext sunctx)
    cdef struct _N_VectorContent_Serial:
        long int length
        realtype *data
        booleantype own_data

    ctypedef _N_VectorContent_Serial *N_VectorContent_Serial


cdef extern from "nvector/nvector_openmp.h":
    cdef N_Vector N_VMake_OpenMP(long int vec_length, realtype *v_data, int num_threads, SUNContext sunctx)

    cdef struct _N_VectorContent_OpenMP:
        long int length
        realtype *data
        booleantype own_data
        int num_threads

    ctypedef _N_VectorContent_OpenMP *N_VectorContent_OpenMP



cdef extern from "cvode/cvode.h":
    int CV_ADAMS
    int CV_BDF
    int CV_FUNCTIONAL
    # int CV_NEWTON
    int CV_NORMAL
    int CV_ONE_STEP

    int CV_SUCCESS
    int CV_TSTOP_RETURN
    int CV_ROOT_RETURN

    int CV_WARNING

    int CV_TOO_MUCH_WORK
    int CV_TOO_MUCH_ACC
    int CV_ERR_FAILURE
    int CV_CONV_FAILURE

    int CV_LINIT_FAIL
    int CV_LSETUP_FAIL
    int CV_LSOLVE_FAIL
    int CV_RHSFUNC_FAIL
    int CV_FIRST_RHSFUNC_ERR
    int CV_REPTD_RHSFUNC_ERR
    int CV_UNREC_RHSFUNC_ERR
    int CV_RTFUNC_FAIL

    int CV_MEM_FAIL
    int CV_MEM_NULL
    int CV_ILL_INPUT
    int CV_NO_MALLOC
    int CV_BAD_K
    int CV_BAD_T
    int CV_BAD_DKY
    int CV_TOO_CLOSE

    ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *user_data)
    ctypedef int (*CVRootFn)(realtype t, N_Vector y, realtype *gout, void *user_data)

    void *CVodeCreate(int lmm, SUNContext sunctx)
    int CVode "CVode"(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, int itask) nogil
    int CVodeSetUserData(void *cvode_mem, void *user_data)
    int CVodeSetMaxOrd(void *cvode_mem, int maxord)
    int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
    int CVodeSetMaxHnilWarns(void *cvode_mem, int mxhnil)
    int CVodeSetStabLimDet(void *cvode_mem, booleantype stldet)
    int CVodeSetInitStep(void *cvode_mem, realtype hin)
    int CVodeSetMinStep(void *cvode_mem, realtype hmin)
    int CVodeSetMaxStep(void *cvode_mem, realtype hmax)
    int CVodeSetStopTime(void *cvode_mem, realtype tstop)
    int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef)
    int CVodeSetMaxNonlinIters(void *cvode_mem, int maxcor)
    int CVodeSetMaxConvFails(void *cvode_mem, int maxncf)
    int CVodeSetNonlinConvCoef(void *cvode_mem, realtype nlscoef)
    int CVodeSetIterType(void *cvode_mem, int iter)
    int CVodeSetRootDirection(void *cvode_mem, int *rootdir)
    int CVodeSetNoInactiveRootWarn(void *cvode_mem)
    int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
    int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0)
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol)
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
    int CVode(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)
    int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky)
    int CVodeGetWorkSpace(void *cvode_mem, long int *lenrw, long int *leniw)
    int CVodeGetNumSteps(void *cvode_mem, long int *nsteps)
    int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevals)
    int CVodeGetNumLinSolvSetups(void *cvode_mem, long int *nlinsetups)
    int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails)
    int CVodeGetLastOrder(void *cvode_mem, int *qlast)
    int CVodeGetCurrentOrder(void *cvode_mem, int *qcur)
    int CVodeGetNumStabLimOrderReds(void *cvode_mem, long int *nslred)
    int CVodeGetActualInitStep(void *cvode_mem, realtype *hinused)
    int CVodeGetLastStep(void *cvode_mem, realtype *hlast)
    int CVodeGetCurrentStep(void *cvode_mem, realtype *hcur)
    int CVodeGetCurrentTime(void *cvode_mem, realtype *tcur)
    int CVodeGetTolScaleFactor(void *cvode_mem, realtype *tolsfac)
    int CVodeGetErrWeights(void *cvode_mem, N_Vector eweight)
    int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele)
    int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals)
    int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
    int CVodeGetIntegratorStats(void *cvode_mem, long int *nsteps,
                                long int *nfevals, long int *nlinsetups,
                                long int *netfails, int *qlast,
                                int *qcur, realtype *hinused, realtype *hlast,
                                realtype *hcur, realtype *tcur)
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters)
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails)
    int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters, long int *nncfails)
    char *CVodeGetReturnFlagName(int flag)
    void CVodeFree(void **cvode_mem)

    int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals)
    int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nrevalsLS)
    int CVodeGetNumJtimesEvals(void *cvode_mem, long int *njevals)

cdef extern from "sunlinsol/sunlinsol_spgmr.h":
    int CVSpgmr(void *cvode_mem, int pretype, int max1)

cdef extern from "cvode/cvode_diag.h":
    int CVDiag(void *cvode_mem)

cdef extern from "cvode/cvode_ls.h":
    # int CVSpilsSetPrecType(void *cvode_mem, int pretype)
    # int CVSpilsSetGSType(void *cvode_mem, int gstype)
    # int CVSpilsSetMaxl(void *cvode_mem, int maxl)
    # int CVSpilsSetEpsLin(void *cvode_mem, realtype eplifac)

    ctypedef int (*CVLsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
                                    booleantype jok, booleantype *jcurPtr,
                                    realtype gamma, void *user_data);
    ctypedef int (*CVLsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
                                    N_Vector r, N_Vector z, realtype gamma,
                                    realtype delta, int lr, void *user_data);
    ctypedef int (*CVLsJacTimesSetupFn)(realtype t, N_Vector y,
                                        N_Vector fy, void *user_data);
    ctypedef int (*CVLsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t,
                                      N_Vector y, N_Vector fy,
                                      void *user_data, N_Vector tmp);

    int CVodeSetPreconditioner(void *cvode_mem, CVLsPrecSetupFn pset, CVLsPrecSolveFn psolve);
    int CVodeSetJacTimes(void *cvode_mem, CVLsJacTimesSetupFn jtsetup, CVLsJacTimesVecFn jtimes);
    # int CVSpilsSetJacTimesVecFn(void *cvode_mem, CVSpilsJacTimesVecFn jtv)

cdef extern from "sundials/sundials_iterative.h":
    int SUN_PREC_NONE
    int SUN_PREC_LEFT
    int SUN_PREC_RIGHT
    int SUN_PREC_BOTH

    int MODIFIED_GS
    int CLASSICAL_GS



cdef extern from "../../atomistic/lib/clib.h":
    void normalise(double * m, int nxyz)


cdef struct cv_userdata:
    void * rhs_fun
    void * y
    void * dm_dt
    void * jvn_fun
    void * v  # is the vector by which the Jacobian must be multiplied
    void * jv  # the value of Jacobian times v

cdef inline int copy_arr2nv(double[:] np_x, N_Vector v):
    cdef long int n = (< N_VectorContent_Serial > v.content).length
    cdef void * data_ptr = &np_x[0]
    memcpy((< N_VectorContent_Serial > v.content).data, data_ptr, n*sizeof(double))
    return 0

cdef inline copy_arr2nv_openmp(double[:] np_x, N_Vector v):
    cdef long int n = (< N_VectorContent_OpenMP > v.content).length
    cdef void * data_ptr = &np_x[0]
    memcpy((< N_VectorContent_OpenMP > v.content).data, data_ptr, n*sizeof(double))
    return 0


cdef inline copy_nv2arr(N_Vector v, double[:] np_x):
    cdef long int n = (< N_VectorContent_Serial > v.content).length
    cdef double * v_data = ( < N_VectorContent_Serial > v.content).data
    memcpy(&np_x[0], v_data, n*sizeof(double))
    return 0

cdef inline copy_nv2arr_openmp(N_Vector v, double[:] np_x):
    cdef long int n = (< N_VectorContent_OpenMP > v.content).length
    cdef double * v_data = ( < N_VectorContent_OpenMP > v.content).data
    memcpy(&np_x[0], v_data, n*sizeof(double))
    return 0


cdef copy_nv2nv(N_Vector v_dest, N_Vector v_src):
    cdef size_t n = (< N_VectorContent_Serial > v_dest.content).length
    cdef double * v_dest_data = ( < N_VectorContent_Serial > v_dest.content).data
    cdef double * v_src_data = ( < N_VectorContent_Serial > v_src.content).data
    memcpy(v_dest_data, v_src_data, n * sizeof(double))
    return 0

cdef copy_nv2nv_openmp(N_Vector v_dest, N_Vector v_src):
    cdef size_t n = (< N_VectorContent_OpenMP > v_dest.content).length
    cdef double * v_dest_data = ( < N_VectorContent_OpenMP > v_dest.content).data
    cdef double * v_src_data = ( < N_VectorContent_OpenMP > v_src.content).data
    memcpy(v_dest_data, v_src_data, n * sizeof(double))
    return 0


cdef int cv_rhs(double t, N_Vector yv, N_Vector yvdot, void * user_data) except -1:
    cdef cv_userdata * ud = <cv_userdata * >user_data
    cdef np.ndarray[double, ndim = 1, mode = 'c'] y_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.y
    cdef np.ndarray ydot_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.dm_dt

    copy_nv2arr(yv, y_arr)
    (< object > ud.rhs_fun)(t, y_arr, ydot_arr)
    copy_arr2nv(ydot_arr, yvdot)
    return 0

cdef int cv_rhs_openmp(double t, N_Vector yv, N_Vector yvdot, void * user_data) except -1:
    cdef cv_userdata * ud = <cv_userdata * >user_data
    cdef np.ndarray[double, ndim = 1, mode = 'c'] y_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.y
    cdef np.ndarray ydot_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.dm_dt

    copy_nv2arr_openmp(yv, y_arr)
    (< object > ud.rhs_fun)(t, y_arr, ydot_arr)
    copy_arr2nv_openmp(ydot_arr, yvdot)
    return 0


cdef int cv_jtimes(N_Vector v, N_Vector Jv, double t, N_Vector y, N_Vector fy, void * user_data, N_Vector tmp) except -1:
    cdef cv_userdata * ud = <cv_userdata * >user_data
    cdef np.ndarray[double, ndim = 1, mode = 'c'] m_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.y
    cdef np.ndarray[double, ndim = 1, mode = 'c'] mp_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.v
    cdef np.ndarray[double, ndim = 1, mode = 'c'] ydot_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.dm_dt
    cdef np.ndarray[double, ndim = 1, mode = 'c'] jv_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.jv

    copy_nv2arr(y, m_arr)
    copy_nv2arr(v, mp_arr)
    copy_nv2arr(fy, ydot_arr)
    (< object > ud.jvn_fun)(mp_arr, jv_arr, t, m_arr, ydot_arr)
    copy_arr2nv(jv_arr, Jv)
    return 0

cdef int cv_jtimes_openmp(N_Vector v, N_Vector Jv, double t, N_Vector y, N_Vector fy, void * user_data, N_Vector tmp) except -1:
    cdef cv_userdata * ud = <cv_userdata * >user_data
    cdef np.ndarray[double, ndim = 1, mode = 'c'] m_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.y
    cdef np.ndarray[double, ndim = 1, mode = 'c'] mp_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.v
    cdef np.ndarray[double, ndim = 1, mode = 'c'] ydot_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.dm_dt
    cdef np.ndarray[double, ndim = 1, mode = 'c'] jv_arr = <np.ndarray[double, ndim = 1, mode = 'c'] > ud.jv

    copy_nv2arr_openmp(y, m_arr)
    copy_nv2arr_openmp(v, mp_arr)
    copy_nv2arr_openmp(fy, ydot_arr)
    (< object > ud.jvn_fun)(mp_arr, jv_arr, t, m_arr, ydot_arr)
    copy_arr2nv_openmp(jv_arr, Jv)
    return 0


# static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
#                   realtype gamma, realtype delta, int lr, void *user_data);
cdef int psolve(double t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                double gamma, double delta, int lr, void * user_data):
    copy_nv2nv(z, r)
    return 0

cdef int psolve_openmp(double t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z,
                       double gamma, double delta, int lr, void * user_data):
    copy_nv2nv_openmp(z, r)
    return 0


cdef class CvodeSolver(object):
    cdef public double t
    cdef public np.ndarray y
    cdef double rtol, atol
    cdef int cvode_already_initialised
    cdef np.ndarray y0
    cdef np.ndarray dm_dt
    cdef np.ndarray mp
    cdef np.ndarray Jmp
    cdef SUNContext sunctx
    cdef N_Vector u_y
    cdef void * cvode_mem
    cdef void * rhs_fun
    cdef void * jvn_fun
    cdef callback_fun
    cdef jtimes_fun
    cdef cv_userdata user_data
    cdef long int nsteps, nfevals, njevals
    cdef int max_num_steps
    cdef int has_jtimes
    cdef str linear_solver
    cdef str parellel_solver
    def __cinit__(self, spins, rhs_fun, jtimes_fun=None, linear_solver="spgmr", rtol=1e-8, atol=1e-8):
        self.t = 0
        self.y0 = spins
        self.dm_dt = np.copy(spins)
        self.y = np.copy(spins)
        self.mp = np.copy(spins)
        self.Jmp = np.copy(spins)

        self.callback_fun = rhs_fun
        self.jtimes_fun = jtimes_fun
        self.jvn_fun = <void * >cv_jtimes
        # wrapper for callback_fun (which is a Python function)
        self.rhs_fun = <void * >cv_rhs
        if linear_solver == "spgmr" or linear_solver == "diag":
            # scaled preconditioned GMRES or diagonal approximate solver
            self.linear_solver = linear_solver
        else:
            raise ValueError(
                "linear_solver is {}, should be spgmr or diag".format(linear_solver))

        self.has_jtimes = 0
        if jtimes_fun is not None:
            self.has_jtimes = 1

        # New in version 6.0.0.
        # All of the SUNDIALS objects (vectors, linear and nonlinear solvers, matrices, etc.)
        # that collectively form a SUNDIALS simulation, hold a reference to a common simulation context object
        # defined by the SUNContext class.
        SUNContext_Create(NULL, & self.sunctx);

        self.user_data = cv_userdata(< void*>self.callback_fun,
                                      < void * >self.y0, < void * >self.dm_dt,
                                      < void * >self.jtimes_fun,
                                      < void * >self.mp, < void * >self.Jmp)

        # The recommended choices for lmm are CV ADAMS for nonstiff problems and CV BDF for
        # stiff problems. The default Newton iteration is recommended for stiff problems, and
        # the fixed-point solver (previously referred to as the functional iteration in this guide) is
        # recommended for nonstiff problems.
        self.cvode_mem = CVodeCreate(CV_BDF, self.sunctx)

        flag = CVodeSetUserData(self.cvode_mem, < void*> & self.user_data)
        self.check_flag(flag, "CVodeSetUserData")

        self.cvode_already_initialised = 0

        self.u_y = N_VMake_Serial(self.y.size, <double *> self.y.data, self.sunctx)

        self.set_initial_value(spins, self.t)
        self.set_options(rtol, atol)

    def reset(self, double[:] spin, t):
        copy_arr2nv(spin, self.u_y)
        CVodeReInit(self.cvode_mem, t, self.u_y)

    def set_initial_value(self, double[:] spin, t):
        self.t = t
        self.y[:] = spin[:]

        cdef double[:] y = self.y
        copy_arr2nv(self.y, self.u_y)

        if self.cvode_already_initialised:
            flag = CVodeReInit(self.cvode_mem, t, self.u_y)
            self.check_flag(flag, "CVodeReInit")
        else:
            flag = CVodeInit(self.cvode_mem, < CVRhsFn > self.rhs_fun, t, self.u_y)
            self.check_flag(flag, "CVodeInit")
            self.cvode_already_initialised = 1

        if self.linear_solver == "diag":
            flag = CVDiag(self.cvode_mem)
            self.check_flag(flag, "CVDiag")
        elif self.linear_solver == "spgmr":
            if self.has_jtimes:
                # The Jacobian preconditioner is set here based on the
                # cvDiurnal_kry.c example from Sundials 6.1.1

                # CVSpgmr(cvode_mem, pretype, maxl) p. 27 of CVODE 2.7 manual
                # flag = CVSpgmr(self.cvode_mem, PREC_LEFT, 300)
                # self.check_flag(flag, "CVSpgmr")

                # Call SUNLinSol_SPGMR to specify the linear solver SPGMR
                # with left preconditioning and the default Krylov dimension
                LS = SUNLinSol_SPGMR(u, SUN_PREC_LEFT, 0, sunctx);
                # TODO:
                # if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

                # Call CVodeSetLinearSolver to attach the linear sovler to CVode
                flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
                self.check_flag(flag, "CVodeSetLinearSolver")
                # if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return 1;

                # Set the Jacobian-times-vector function */
                flag = CVodeSetJacTimes(cvode_mem, NULL, jtv);
                self.check_flag(flag, "CVodeSetJacTimes")

                # functions below in p. 37 CVODE 2.7 manual
                # flag = CVodeSetJacTimes(self.cvode_mem, NULL, < CVSpilsJacTimesVecFn > self.jvn_fun)
                # self.check_flag(flag, "CVSpilsSetJacTimesVecFn")

                flag = CVodeSetPreconditioner(self.cvode_mem,
                                             < CVLsPrecSetupFn > self.Precond,
                                             < CVLsPrecSolveFn > psolve)
                self.check_flag(flag, "CVodeSetPreconditioner")
            else:
                # this will use the SPGMR without preconditioner and without
                # our computation of the product J * m'. Instead, it uses
                # a difference quotient approximation of the product.
                # c.f. Sec 4.6.7 in CVODE manual
                # Actually, it's the same Jacobian approximation as used
                # in CVDiag (only difference is CVDiag is a direct linear
                # solver).
                flag = SUNLinSol_SPGMR(self.cvode_mem, SUN_PREC_NONE, 300)
                # TODO:
                # self.check_flag(flag, "SUNLinSol_SPGMR")
        else:
            raise RuntimeError(
                "linear_solver is {}, should be spgmr or diag".format(self.linear_solver))

    def set_options(self, rtol, atol, max_num_steps=100000, max_ord=None):
        self.rtol = rtol
        self.atol = atol
        self.max_num_steps = max_num_steps

        # Set tolerances
        flag = CVodeSStolerances(self.cvode_mem, self.rtol, self.atol)

        # maximum order of the linear multistep method
        # c.f. Suess et al. / Journal of Magnetism and Magnetic Materials 248 (2002) 298–311
        # where it speeds up time integration considerably (Fig. 3)
        # default is 5 for BDF method, 2 is a good alternative value to try
        if max_ord is not None:
            flag = CVodeSetMaxOrd(self.cvode_mem, max_ord)
            self.check_flag(flag, "CVodeSetMaxOrd")

        # Set maximum number of iteration steps (?)
        flag = CVodeSetMaxNumSteps(self.cvode_mem, max_num_steps)
        self.check_flag(flag, "CVodeSetMaxNumSteps")

        flag = CVodeReInit(self.cvode_mem, self.t, self.u_y)
        self.check_flag(flag, "CVodeReInit")

    # TODO: Instead of using flags we should use Sundials internal flag
    # like CV_SUCCESS
    cpdef int run_until(self, double t_final) except -1:
        cdef int flag
        cdef double t_returned
        flag = CVode(self.cvode_mem, t_final, self.u_y, & t_returned, CV_NORMAL)
        self.check_flag(flag, "CVode")
        self.t = t_returned
        return 0

    # From exmaples: cvDiurnal_kry.c in Sundials repo:
    # static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
    #                    booleantype *jcurPtr, realtype gamma, void *user_data)
    cdef int Precond(self, double t, N_Vector y, N_Vector fy, booleantype jok,
                     booleantype * jcurPtr, double gamma, void * user_data):
        if not jok:
            copy_nv2arr(y, self.y)
        return 0

    def check_flag(self, flag, fun_name):
        if flag != 0:
            raise RuntimeError("CVODE function {} failed!".format(fun_name))

    def rhs_evals(self):
        CVodeGetNumRhsEvals(self.cvode_mem, & self.nfevals)
        return self.nfevals

    def stat(self):
        CVodeGetNumSteps(self.cvode_mem, & self.nsteps)
        CVodeGetNumRhsEvals(self.cvode_mem, & self.nfevals)
        CVodeGetNumJtimesEvals(self.cvode_mem, & self.njevals)
        return self.nsteps, self.nfevals, self.njevals

    def get_current_step(self):
        cdef double step
        CVodeGetCurrentStep(self.cvode_mem, & step)
        return step

    def __repr__(self):
        return "nsteps = {}, nfevals = {}, njevals = {}".format(
            self.nsteps, self.nfevals, self.njevals)

    def __str__(self):
        return "{} ({})".format(self.__class__.__name__,
                                self.__repr__())

    def __dealloc__(self):
        self.user_data.rhs_fun = NULL
        self.user_data.y = NULL
        self.user_data.jvn_fun = NULL
        self.user_data.dm_dt = NULL
        self.user_data.v = NULL
        self.user_data.jv = NULL
        N_VDestroy_Serial(self.u_y)
        CVodeFree(& self.cvode_mem)
        SUNContext_Free(& self.sunctx)

cdef class CvodeSolver_OpenMP(object):
    cdef public double t
    cdef public np.ndarray y
    cdef double rtol, atol
    cdef int cvode_already_initialised
    cdef np.ndarray y0
    cdef np.ndarray dm_dt
    cdef np.ndarray mp
    cdef np.ndarray Jmp
    cdef SUNContext sunctx
    cdef N_Vector u_y
    cdef void * cvode_mem
    cdef void * rhs_fun
    cdef void * jvn_fun
    cdef callback_fun
    cdef jtimes_fun
    cdef cv_userdata user_data
    cdef long int nsteps, nfevals, njevals
    cdef int max_num_steps
    cdef int has_jtimes
    cdef str linear_solver
    cdef str parellel_solver
    cdef int num_threads

    def __cinit__(self, spins, rhs_fun, jtimes_fun=None, linear_solver="spgmr", rtol=1e-8, atol=1e-8):
        self.num_threads = openmp.omp_get_max_threads()
        print("Number of threads (CVODE) = {}".format(self.num_threads))
        self.t = 0
        self.y0 = spins
        self.dm_dt = np.copy(spins)
        self.y = np.copy(spins)
        self.mp = np.copy(spins)
        self.Jmp = np.copy(spins)

        self.callback_fun = rhs_fun
        self.jtimes_fun = jtimes_fun
        self.jvn_fun = <void * >cv_jtimes_openmp
        # wrapper for callback_fun (which is a Python function)
        self.rhs_fun = <void * >cv_rhs_openmp
        if linear_solver == "spgmr" or linear_solver == "diag":
            # scaled preconditioned GMRES or diagonal approximate solver
            self.linear_solver = linear_solver
        else:
            raise ValueError(
                "linear_solver is {}, should be spgmr or diag".format(linear_solver))

        self.has_jtimes = 0
        if jtimes_fun is not None:
            self.has_jtimes = 1

        self.user_data = cv_userdata(< void*>self.callback_fun,
                                      < void * >self.y0, < void * >self.dm_dt,
                                      < void * >self.jtimes_fun,
                                      < void * >self.mp, < void * >self.Jmp)

        SUNContext_Create(NULL, & self.sunctx);
        # Newton iterator is set by default now (Sundials 4.0)
        self.cvode_mem = CVodeCreate(CV_BDF, self.sunctx)

        flag = CVodeSetUserData(self.cvode_mem, < void*> & self.user_data)
        self.check_flag(flag, "CVodeSetUserData")

        self.cvode_already_initialised = 0

        self.u_y = N_VMake_OpenMP(self.y.size, <double *> self.y.data, self.num_threads, self.sunctx)

        self.set_initial_value(spins, self.t)
        self.set_options(rtol, atol)

    def reset(self, np.ndarray[double, ndim=1, mode="c"] spin, t):
        copy_arr2nv_openmp(spin, self.u_y)
        CVodeReInit(self.cvode_mem, t, self.u_y)

    def set_initial_value(self, np.ndarray[double, ndim=1, mode="c"] spin, t):
        self.t = t
        self.y[:] = spin[:]

        cdef np.ndarray[double, ndim = 1, mode = "c"] y = self.y
        copy_arr2nv_openmp(self.y, self.u_y)

        if self.cvode_already_initialised:
            flag = CVodeReInit(self.cvode_mem, t, self.u_y)
            self.check_flag(flag, "CVodeReInit")
        else:
            flag = CVodeInit(self.cvode_mem, < CVRhsFn > self.rhs_fun, t, self.u_y)
            self.check_flag(flag, "CVodeInit")
            self.cvode_already_initialised = 1

        if self.linear_solver == "diag":
            flag = CVDiag(self.cvode_mem)
            self.check_flag(flag, "CVDiag")
        elif self.linear_solver == "spgmr":
            if self.has_jtimes:

                LS = SUNLinSol_SPGMR(u, SUN_PREC_LEFT, 0, sunctx);

                flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
                self.check_flag(flag, "CVodeSetLinearSolver")

                flag = CVodeSetJacTimes(cvode_mem, NULL, jtv);
                self.check_flag(flag, "CVodeSetJacTimes")

                flag = CVodeSetPreconditioner(self.cvode_mem,
                                             < CVLsPrecSetupFn > self.Precond,
                                             < CVLsPrecSolveFn > psolve_openmp)
                self.check_flag(flag, "CVodeSetPreconditioner")

            else:
                # this will use the SPGMR without preconditioner and without
                # our computation of the product J * m'. Instead, it uses
                # a difference quotient approximation of the product.
                # c.f. Sec 4.6.7 in CVODE manual
                # Actually, it's the same Jacobian approximation as used
                # in CVDiag (only difference is CVDiag is a direct linear
                # solver).
                flag = SUNLinSol_SPGMR(self.cvode_mem, SUN_PREC_NONE, 300)
                # self.check_flag(flag, "CVSpgmr")
        else:
            raise RuntimeError(
                "linear_solver is {}, should be spgmr or diag".format(self.linear_solver))

    def set_options(self, rtol, atol, max_num_steps=100000, max_ord=None):
        self.rtol = rtol
        self.atol = atol
        self.max_num_steps = max_num_steps

        # Set tolerances
        flag = CVodeSStolerances(self.cvode_mem, self.rtol, self.atol)

        # maximum order of the linear multistep method
        # c.f. Suess et al. / Journal of Magnetism and Magnetic Materials 248 (2002) 298–311
        # where it speeds up time integration considerably (Fig. 3)
        # default is 5 for BDF method, 2 is a good alternative value to try
        if max_ord is not None:
            flag = CVodeSetMaxOrd(self.cvode_mem, max_ord)
            self.check_flag(flag, "CVodeSetMaxOrd")

        # Set maximum number of iteration steps (?)
        flag = CVodeSetMaxNumSteps(self.cvode_mem, max_num_steps)
        self.check_flag(flag, "CVodeSetMaxNumSteps")

        flag = CVodeReInit(self.cvode_mem, self.t, self.u_y)
        self.check_flag(flag, "CVodeReInit")

    cpdef int run_until(self, double t_final) except -1:
        cdef int flag
        cdef double t_returned
        flag = CVode(self.cvode_mem, t_final, self.u_y, & t_returned, CV_NORMAL)
        self.check_flag(flag, "CVode")
        self.t = t_returned
        return 0

    cdef int Precond(self, double t, N_Vector y, N_Vector fy,
                     booleantype jok, booleantype * jcurPtr, double gamma,
                     void * user_data):
        if not jok:
            copy_nv2arr_openmp(y, self.y)
        return 0

    def check_flag(self, flag, fun_name):
        if flag != 0:
            raise RuntimeError("CVODE function {} failed!".format(fun_name))

    def rhs_evals(self):
        CVodeGetNumRhsEvals(self.cvode_mem, & self.nfevals)
        return self.nfevals

    def stat(self):
        CVodeGetNumSteps(self.cvode_mem, & self.nsteps)
        CVodeGetNumRhsEvals(self.cvode_mem, & self.nfevals)
        CVodeGetNumJtimesEvals(self.cvode_mem, & self.njevals)
        return self.nsteps, self.nfevals, self.njevals

    def get_current_step(self):
        cdef double step
        CVodeGetCurrentStep(self.cvode_mem, & step)
        return step

    def __repr__(self):
        return "nsteps = {}, nfevals = {}, njevals = {}".format(
            self.nsteps, self.nfevals, self.njevals)

    def __str__(self):
        return "{} ({})".format(self.__class__.__name__,
                                self.__repr__())

    def __dealloc__(self):
        self.user_data.rhs_fun = NULL
        self.user_data.y = NULL
        self.user_data.jvn_fun = NULL
        self.user_data.dm_dt = NULL
        self.user_data.v = NULL
        self.user_data.jv = NULL
        N_VDestroy_OpenMP(self.u_y)
        CVodeFree(& self.cvode_mem)
        SUNContext_Free(& self.sunctx)
