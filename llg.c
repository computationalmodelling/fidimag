#include "clib.h"

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

		mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
		relax = c * (1 - mm);
		dm_dt[i] += relax * m[i];
		dm_dt[j] += relax * m[j];
		dm_dt[k] += relax * m[k];

	}

}

