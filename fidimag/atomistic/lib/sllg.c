#include "clib.h"
#include "llg_random.h"


ode_solver *create_ode_plan(void) {

	ode_solver *plan = (ode_solver*) malloc(sizeof(ode_solver));

	return plan;
}


void init_solver(ode_solver *s, double k_B, double theta, int nxyz, double dt, double gamma) {

	s->theta = theta;
	s->theta1 = 1-0.5/theta;
	s->theta2 = 0.5/theta;

	s->dt = dt;
	s->nxyz = nxyz;
	s->gamma = gamma;
	
	s->Q = 2 * k_B * dt / gamma;

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


void llg_rhs_dw(ode_solver *s, double *m, double *h, double *dm, double *T, double *alpha, double *mu_s_inv, int *pins) {

	int i, j, k;

	double mth0, mth1, mth2;
	
	int nxyz = s->nxyz;
	double *eta = &s->eta[0];
	double dt = s->dt;
	double Q=s->Q;

	double coeff, q;
	double hi,hj,hk;

	#pragma omp parallel for private(i,j,k,q,coeff,hi,hj,hk,mth0,mth1,mth2)
	for (int id = 0; id < nxyz; id++) {
		i = 3*id;
		j = i+1;
		k = j+1;
		
		if (pins[i]>0){
			 dm[i] = 0;
			 dm[j] = 0;
			 dm[k] = 0;
			 continue;
		}


        coeff = -s->gamma/ (1.0 + alpha[id] * alpha[id]) ;
		q = sqrt(Q * alpha[i] * T[i] * mu_s_inv[id]);
		
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

void run_step1(ode_solver *s, double *m, double *h, double *m_pred, double *T, double *alpha, double *mu_s_inv, int *pins) {
	int i;
	int nxyz = s->nxyz;
	double *dm1 = s->dm1;
	double theta = s->theta;

	gauss_random_vec(s->eta, 3 * nxyz);
	llg_rhs_dw(s, m, h, dm1, T, alpha, mu_s_inv, pins);

	//#pragma omp parallel for private(i)
	for (i = 0; i < 3 * s->nxyz; i++) {
		m_pred[i] = m[i] + theta * dm1[i];
	}

}

void run_step2(ode_solver *s, double *m_pred, double *h, double *m, double *T, double *alpha, double *mu_s_inv, int *pins) {
	int i, j, k;
	int nxyz = s->nxyz;
	double *dm1 = s->dm1;
	double *dm2 = s->dm2;
	double theta1 = s->theta1;
	double theta2 = s->theta2;

	llg_rhs_dw(s, m_pred, h, dm2, T, alpha, mu_s_inv, pins);

	//#pragma omp parallel for private(i)
	for (i = 0; i < 3 * nxyz; i++) {
		m[i] += (theta1 * dm1[i] + theta2 * dm2[i]);
	}

	double mm;
	#pragma omp parallel for private(i,j,k,mm)
	for (int id = 0; id < s->nxyz; id++) {

		if (pins[id]>0){
			continue;
		}

		i= 3*id;
		j = i + 1;
		k = j + 1;

		mm = 1.0 / sqrt(m[i] * m[i] + m[j] * m[j] + m[k] * m[k]);
		m[i] *= mm;
		m[j] *= mm;
		m[k] *= mm;
	}
	//we have to say that this kind of method is quite inaccurate, so in future we can try other methods

}

void finalize_ode_plan(ode_solver *plan) {
	free(plan->dm1);
	free(plan->dm2);
	free(plan->eta);
	free(plan);
}


void normalise(double *m, int nxyz){
	int i, j, k;
	double mm;
	for (int id = 0; id < nxyz; id++) {
			i = 3*id;
			j = i + 1;
			k = j + 1;

			mm = 1.0 / sqrt(m[i] * m[i] + m[j] * m[j] + m[k] * m[k]);
			m[i] *= mm;
			m[j] *= mm;
			m[k] *= mm;

	}
}
