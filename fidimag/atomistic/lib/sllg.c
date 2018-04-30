#include "clib.h"

/*
* n is the spin number
* eta is the random number array
*/
void llg_rhs_dw_c(double *restrict m, double *restrict h, double *restrict dm, double *restrict T, double *restrict alpha,
                  double *restrict mu_s_inv, int *restrict pins, double *restrict eta, int n, double gamma, double dt) {
        
        double k_B = 1.3806505e-23;
        double Q = 2 * k_B * dt / gamma;

	//#pragma omp parallel for
	for (int id = 0; id < n; id++) {
		int i = 3*id;
		int j = i+1;
		int k = j+1;
		
		if (pins[id]>0){
			 dm[i] = 0;
			 dm[j] = 0;
			 dm[k] = 0;
			 continue;
		}


        double coeff = -gamma/ (1.0 + alpha[id] * alpha[id]);
		double q = sqrt(Q * alpha[id] * T[id] * mu_s_inv[id]);
		
		double hi = h[i]*dt + eta[i]*q;
		double hj = h[j]*dt + eta[j]*q;
		double hk = h[k]*dt + eta[k]*q;
		
		double mth0 = coeff * (m[j] * hk - m[k] * hj);
		double mth1 = coeff * (m[k] * hi - m[i] * hk);
		double mth2 = coeff * (m[i] * hj - m[j] * hi);

		dm[i] = mth0 + alpha[id] * (m[j] * mth2 - m[k] * mth1);
		dm[j] = mth1 + alpha[id] * (m[k] * mth0 - m[i] * mth2);
		dm[k] = mth2 + alpha[id] * (m[i] * mth1 - m[j] * mth0);
		
	}
}

void normalise(double *restrict m, int *restrict pins, int n){
	int i, j, k;
	double mm;
	for (int id = 0; id < n; id++) {
			i = 3*id;
			j = i + 1;
			k = j + 1;

		        if (pins[id]>0) continue;
	

			mm = 1.0 / sqrt(m[i] * m[i] + m[j] * m[j] + m[k] * m[k]);
			m[i] *= mm;
			m[j] *= mm;
			m[k] *= mm;

	}
}
