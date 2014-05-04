#include "clib.h"
#include "llg_random.h"


ode_solver *create_ode_plan() {

	ode_solver *plan = (ode_solver*) malloc(sizeof(ode_solver));

	return plan;
}

//used for debug
void print_g(char *str, double *x, int n) {
	int i;
	printf("%s:\n         ", str);
	for (i = 0; i < n; i++) {
		printf("%g ", x[i]);
	}
	printf("\n");

}

void llg_rhs_dw(ode_solver *s, double *m, double *h, double *dm, double *T, double *alpha) {

	int i, j, k;

	double mth0, mth1, mth2;
	
	int nxyz = s->nxyz;
	double *eta = &s->eta[0];
	double dt = s->dt;

	double coeff = s->coeff;
	double q, alpha_inv;
	double hi,hj,hk;


	for (i = 0; i < nxyz; i++) {

		j = i + nxyz;
		k = j + nxyz;
		
		alpha_inv = 1.0/ (1.0 + alpha[i] * alpha[i]);
        coeff = -s->gamma * alpha_inv ;
		q = sqrt(s->Q * alpha[i] *alpha_inv * T[i]);
		
		hi = h[i]*dt + eta[i]*q;
		hj = h[j]*dt + eta[j]*q;
		hk = h[k]*dt + eta[k]*q;
		
		mth0 = coeff * (m[j] * hk - m[k] * hj);
		mth1 = coeff * (m[k] * hi - m[i] * hk);
		mth2 = coeff * (m[i] * hj - m[j] * hi);

		dm[i] = mth0 + alpha[i] * (m[j] * mth2 - m[k] * mth1);
		dm[j] = mth1 + alpha[i] * (m[k] * mth0 - m[i] * mth2);
		dm[k] = mth2 + alpha[i] * (m[i] * mth1 - m[j] * mth0);
		
	}
}

void init_solver(ode_solver *s, double mu_s, int nxyz, double dt, double gamma) {

	s->theta = 2.0 / 3.0;
	s->theta1 = 0.25;
	s->theta2 = 0.75;

	s->dt = dt;
	s->nxyz = nxyz;
	s->gamma = gamma;
	
	double k_B = 1.3806505e-23;
	s->Q = 2 * k_B / (gamma * mu_s) * dt;

	s->dm1 = (double*) malloc(3 * nxyz * sizeof(double));
	s->dm2 = (double*) malloc(3 * nxyz * sizeof(double));
	s->eta = (double*) malloc(3 * nxyz * sizeof(double));

	int i = 0;
	for (i = 0; i < 3 * nxyz; i++) {
		s->dm1[i] = 0;
		s->dm2[i] = 0;
		s->eta[i] = 0;
	}

	initial_random();
}

void run_step1(ode_solver *s, double *m, double *h, double *m_pred, double *T, double *alpha) {
	int i, j, k;
	int nxyz = s->nxyz;
	double *dm1 = s->dm1;
	double theta = s->theta;

	gauss_random_vec(s->eta, 3 * nxyz);
	llg_rhs_dw(s, m, h, dm1, T, alpha);

	for (i = 0; i < 3 * s->nxyz; i++) {
		m_pred[i] = m[i] + theta * dm1[i];
	}

}

void run_step2(ode_solver *s, double *m_pred, double *h, double *m, double *T, double *alpha) {
	int i, j, k;
	int nxyz = s->nxyz;
	double *dm1 = s->dm1;
	double *dm2 = s->dm2;
	double theta1 = s->theta1;
	double theta2 = s->theta2;

	llg_rhs_dw(s, m_pred, h, dm2, T, alpha);

	for (i = 0; i < 3 * nxyz; i++) {
		m[i] += (theta1 * dm1[i] + theta2 * dm2[i]);
	}

	double mm;
	for (i = 0; i < s->nxyz; i++) {
		j = i + nxyz;
		k = j + nxyz;
		mm = 1.0 / sqrt(m[i] * m[i] + m[j] * m[j] + m[k] * m[k]);
		m[i] *= mm;
		m[j] *= mm;
		m[k] *= mm;
	}

}

void finalize_ode_plan(ode_solver *plan) {
	free(plan->dm1);
	free(plan->dm2);
	free(plan->eta);
	free(plan);
}


int check_array(double *a, double *b, int n){
    int i;
    for(i=0;i<n;i++){
        if (fabs(a[i]-b[i])>0){
            printf("a=%g   b=%g\n",a[i],b[i]);
            return -1;
        }
    }
    return 0;

}

void llg_rhs(double *dm_dt, double *m, double *h, double *alpha, double gamma, int nxyz) {

	int i, j, k;

	double mth0, mth1, mth2;
    double coeff, mm, mh, c;
    double hpi,hpj,hpk;

	for (i = 0; i < nxyz; i++) {
		j = i + nxyz;
		k = j + nxyz;

		coeff = -gamma/(1+alpha[i]*alpha[i]);
        
        mm = m[i]*m[i] + m[j]*m[j] + m[k]*m[k];
        mh = m[i]*h[i] + m[j]*h[j] + m[k]*h[k];
        
        //suppose m is normalised, i.e., mm=1; hp=mm.h-mh.m=-mx(mxh)
        //hpi = mm*h[i] - mh*m[i];
        //hpj = mm*h[j] - mh*m[j];
        //hpk = mm*h[k] - mh*m[k];
        hpi = h[i] - mh*m[i];
        hpj = h[j] - mh*m[j];
        hpk = h[k] - mh*m[k];
        
        mth0 = (m[j] * hpk - m[k] * hpj);
		mth1 = (m[k] * hpi - m[i] * hpk);
		mth2 = (m[i] * hpj - m[j] * hpi);
        
        dm_dt[i] = coeff*(mth0 - hpi * alpha[i]);
        dm_dt[j] = coeff*(mth1 - hpj * alpha[i]);
        dm_dt[k] = coeff*(mth2 - hpk * alpha[i]);
        
        // in future, we will try the new method to integrate the LLG equation,
        // A mixed mid-point Runge-Kutta like scheme for the integration of Landau-Lifshitz equation
        // Journal of Applied Physics 115, 17D101 (2014)
        // if possible, we can combine it with adaptive step size, don't know how to do but it's worth a try.
        c = 6*sqrt(dm_dt[i]*dm_dt[i]+dm_dt[j]*dm_dt[j]+dm_dt[k]*dm_dt[k]);
        dm_dt[i] += c*(1-mm)*m[i];
        dm_dt[j] += c*(1-mm)*m[j];
        dm_dt[k] += c*(1-mm)*m[k];

	}
    
}


void llg_s_rhs(double *dm_dt, double *m, double *h, double *alpha, double *chi, double gamma, int nxyz) {
    
	int i, j, k;
    
	double mth0, mth1, mth2;
	double coeff = - gamma;
    double mm, c;
    double hpi,hpj,hpk,mh;
    
	for (i = 0; i < nxyz; i++) {
		j = i + nxyz;
		k = j + nxyz;
        
        if (chi[i] > 0.0){
            mth0 = coeff * (m[j] * h[k] - m[k] * h[j]);
            mth1 = coeff * (m[k] * h[i] - m[i] * h[k]);
            mth2 = coeff * (m[i] * h[j] - m[j] * h[i]);
        
            dm_dt[i] = mth0 + alpha[i] * gamma * h[i];
            dm_dt[j] = mth1 + alpha[i] * gamma * h[j];
            dm_dt[k] = mth2 + alpha[i] * gamma * h[k];
        
            mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
            c = (1-mm)/chi[i];
        
            dm_dt[i] += c*m[i];
            dm_dt[j] += c*m[j];
            dm_dt[k] += c*m[k];
        
        }else{ //do LL equation, simliar to the function of llg_rhs
   
            mh = m[i]*h[i] + m[j]*h[j] + m[k]*h[k];
            
            hpi = h[i] - mh*m[i];
            hpj = h[j] - mh*m[j];
            hpk = h[k] - mh*m[k];
            
            mth0 = (m[j] * hpk - m[k] * hpj);
            mth1 = (m[k] * hpi - m[i] * hpk);
            mth2 = (m[i] * hpj - m[j] * hpi);
            
            dm_dt[i] = coeff*(mth0 - hpi * alpha[i]);
            dm_dt[j] = coeff*(mth1 - hpj * alpha[i]);
            dm_dt[k] = coeff*(mth2 - hpk * alpha[i]);
        }
         
	}
    
}


void normalise(double *m, int nxyz){
	int i, j, k;
	double mm;
	for (i = 0; i < nxyz; i++) {
			j = i + nxyz;
			k = j + nxyz;

			mm = 1.0 / sqrt(m[i] * m[i] + m[j] * m[j] + m[k] * m[k]);
			m[i] *= mm;
			m[j] *= mm;
			m[k] *= mm;

	}
}
