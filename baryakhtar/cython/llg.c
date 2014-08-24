#include "baryakhtar_clib.h"

void llg_rhs_baryakhtar(double *dm_dt, double *m, double *h, double *delta_h,
		double *alpha, double beta, int *pins,
		double gamma, int nxyz, int do_procession) {

	#pragma omp parallel for
	for (int i = 0; i < nxyz; i++) {
		int j = i + nxyz;
		int k = j + nxyz;

		if (pins[i]>0){
			 dm_dt[i] = 0;
			 dm_dt[j] = 0;
			 dm_dt[k] = 0;
			 continue;
		}

		double coeff = -gamma;

        if (do_procession){
        	dm_dt[i] = coeff*cross_x(m[i],m[j],m[k],h[i],h[j],h[k]);
        	dm_dt[j] = coeff*cross_y(m[i],m[j],m[k],h[i],h[j],h[k]);
        	dm_dt[k] = coeff*cross_z(m[i],m[j],m[k],h[i],h[j],h[k]);
        }
        

    	dm_dt[i] += gamma*(alpha[i]*h[i] - beta*delta_h[i]);
    	dm_dt[j] += gamma*(alpha[i]*h[j] - beta*delta_h[j]);
    	dm_dt[k] += gamma*(alpha[i]*h[k] - beta*delta_h[k]);

	}
    
}


