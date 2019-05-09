#include "baryakhtar_clib.h"

void llg_rhs_baryakhtar(double *dm_dt, double *m, double *h, double *delta_h,
		double *alpha, double beta, int *pins,
		double gamma, int nxyz, int do_precession) {

	#pragma omp parallel for
	for (int id = 0; id < nxyz; id++) {
        int i = 3*id;
		int j = i + 1;
		int k = j + 1;

		if (pins[id]>0){
			 dm_dt[i] = 0;
			 dm_dt[j] = 0;
			 dm_dt[k] = 0;
			 continue;
		}

		double coeff = -gamma;

        if (do_precession){
        	dm_dt[i] = coeff*cross_x(m[i],m[j],m[k],h[i],h[j],h[k]);
        	dm_dt[j] = coeff*cross_y(m[i],m[j],m[k],h[i],h[j],h[k]);
        	dm_dt[k] = coeff*cross_z(m[i],m[j],m[k],h[i],h[j],h[k]);
        }

    	dm_dt[i] += gamma*(alpha[i]*h[i] - beta*delta_h[i]);
    	dm_dt[j] += gamma*(alpha[i]*h[j] - beta*delta_h[j]);
    	dm_dt[k] += gamma*(alpha[i]*h[k] - beta*delta_h[k]);

	}

}

void llg_rhs_baryakhtar_reduced(double *dm_dt, double *m, double *hp, double *delta_hp,
                        double *alpha, double beta, int *pins,
                        double gamma, int nxyz, int do_precession, double default_c) {

    #pragma omp parallel for
	for (int id = 0; id < nxyz; id++) {
        int i = 3*id;
		int j = i + 1;
		int k = j + 1;

		if (pins[id]>0){
            dm_dt[i] = 0;
            dm_dt[j] = 0;
            dm_dt[k] = 0;
            continue;
		}

		double coeff = -gamma;

        if (do_precession){
        	dm_dt[i] = coeff*cross_x(m[i],m[j],m[k],hp[i],hp[j],hp[k]);
        	dm_dt[j] = coeff*cross_y(m[i],m[j],m[k],hp[i],hp[j],hp[k]);
        	dm_dt[k] = coeff*cross_z(m[i],m[j],m[k],hp[i],hp[j],hp[k]);
        }

        double hpx = alpha[i]*hp[i]-beta*delta_hp[i];
        double hpy = alpha[i]*hp[j]-beta*delta_hp[j];
        double hpz = alpha[i]*hp[k]-beta*delta_hp[k];

        double mm = m[i]*m[i] + m[j]*m[j] + m[k]*m[k];
        double mh = m[i]*hpx + m[j]*hpy + m[k]*hpz;

        //suppose m is normalised, i.e., hp=mm.h-mh.m=-mx(mxh)

    	dm_dt[i] += gamma*(mm*hpx - mh*m[i]);
    	dm_dt[j] += gamma*(mm*hpy - mh*m[j]);
    	dm_dt[k] += gamma*(mm*hpz - mh*m[k]);

        double c=0;
        if (default_c<0){
        	c = 6*sqrt(dm_dt[i]*dm_dt[i]+dm_dt[j]*dm_dt[j]+dm_dt[k]*dm_dt[k]);
        }else{
        	c = default_c;
        }
        //printf("%0.15g   %0.15g\n", c, default_c);

        dm_dt[i] += c*(1-mm)*m[i];
        dm_dt[j] += c*(1-mm)*m[j];
        dm_dt[k] += c*(1-mm)*m[k];

	}

}
