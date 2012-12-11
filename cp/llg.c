#include "clib.h"
#include "llg_random.h"

void llg_rhs(double *dm_dt, double *m, double *h, double gamma, double alpha,
		double mu_s, int nxyz, double c) {

	int i, j, k;

	double mth0, mth1, mth2;
	double coeff = -gamma / (1 + alpha * alpha) / mu_s;
	double mm, relax;

	for (i = 0; i < nxyz; i++) {
		j = i + nxyz;
		k = j + nxyz;

		mth0 = coeff * (m[j] * h[k] - m[k] * h[j]);
		mth1 = coeff * (m[k] * h[i] - m[i] * h[k]);
		mth2 = coeff * (m[i] * h[j] - m[j] * h[i]);

		dm_dt[i] = mth0 + alpha * (m[j] * mth2 - m[k] * mth1);
		dm_dt[j] = mth1 + alpha * (m[k] * mth0 - m[i] * mth2);
		dm_dt[k] = mth2 + alpha * (m[i] * mth1 - m[j] * mth0);

		mm = 1.0 / sqrt(m[i] * m[i] + m[j] * m[j] + m[k] * m[k]);
		m[i] *= mm;
		m[j] *= mm;
		m[k] *= mm;

		/*
		 relax = c * (1 - mm);
		 dm_dt[i] += relax * m[i];
		 dm_dt[j] += relax * m[j];
		 dm_dt[k] += relax * m[k];
		*/
	}

}

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

void llg_rhs_dw(ode_solver *s, double *m, double *h, double *dm) {

	int i, j, k;

	double mth0, mth1, mth2;
	//double coeff = -gamma / (1 + alpha * alpha) / mu_s;
	double mm, relax;
	double sqrt_dt = sqrt(s->dt);

	int nxyz = s->nxyz;
	double *eta = s->eta;
	double dt = s->dt;
	double Q = s->Q;
	double coeff = s->coeff;
	double alpha = s->alpha;

	gauss_random_vec(eta, 3 * s->nxyz, sqrt_dt);

	for (i = 0; i < nxyz; i++) {
		j = i + nxyz;
		k = j + nxyz;

		mth0 = coeff * (m[j] * h[k] - m[k] * h[j]) * dt;
		mth1 = coeff * (m[k] * h[i] - m[i] * h[k]) * dt;
		mth2 = coeff * (m[i] * h[j] - m[j] * h[i]) * dt;

		mth0 += coeff * (m[j] * eta[k] - m[k] * eta[j]) * Q;
		mth1 += coeff * (m[k] * eta[i] - m[i] * eta[k]) * Q;
		mth2 += coeff * (m[i] * eta[j] - m[j] * eta[i]) * Q;

		dm[i] = mth0 + alpha * (m[j] * mth2 - m[k] * mth1);
		dm[j] = mth1 + alpha * (m[k] * mth0 - m[i] * mth2);
		dm[k] = mth2 + alpha * (m[i] * mth1 - m[j] * mth0);

		/*
		 mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
		 relax = s->c * (1 - mm);
		 dm[i] += relax * m[i] * dt;
		 dm[j] += relax * m[j] * dt;
		 dm[k] += relax * m[k] * dt;
		 */

	}
}

void init_solver(ode_solver *s, double mu_s, int nxyz, double dt, double gamma,
		double alpha, double T, double c) {

	s->theta = 2.0 / 3.0;
	s->theta1 = 1.0 - 0.5 / s->theta;
	s->theta2 = 0.5 / s->theta;

	s->T = T;
	s->dt = dt;
	s->nxyz = nxyz;
	s->alpha = alpha;
	s->gamma = gamma;
	s->coeff = -gamma / (1.0 + alpha * alpha) / 1.0;

	double k_B = 1.3806505e-23;
	//double mu_0 = 4*M_PI*1e-7;
	s->Q = sqrt(2 * k_B * alpha * T / (gamma * mu_s));
	s->c = c;

	/*
	 printf("nxyz=%d dt=%g coeff=%g  Q=%g  Q*sqrt(dt)=%g\n", nxyz, dt, s->coeff,
	 s->Q, s->Q * sqrt(dt));
	 */

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

void run_step1(ode_solver *s, double *m, double *h, double *m_pred) {
	int i;
	double *dm1 = s->dm1;
	double theta = s->theta;

	llg_rhs_dw(s, m, h, dm1);

	for (i = 0; i < 3 * s->nxyz; i++) {
		m_pred[i] = m[i] + theta * dm1[i];
	}

}

void run_step2(ode_solver *s, double *m_pred, double *h, double *m) {
	int i, j, k;
	int nxyz = s->nxyz;
	double *dm1 = s->dm1;
	double *dm2 = s->dm2;
	double theta1 = s->theta1;
	double theta2 = s->theta2;

	llg_rhs_dw(s, m_pred, h, dm2);

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
