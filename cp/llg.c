#include "clib.h"

void llg_rhs(double *dm_dt, double *m, double *h, double *alpha, int *pins, double gamma, int nxyz) {

	int i, j, k;

	double mth0, mth1, mth2;
    double coeff, mm, mh, c;
    double hpi,hpj,hpk;

	#pragma omp parallel for private(i,j,k,coeff,mm, mh, c, mth0, mth1, mth2, hpi,hpj,hpk)
	for (i = 0; i < nxyz; i++) {
		j = i + nxyz;
		k = j + nxyz;

		if (pins[i]>0){
			 dm_dt[i] = 0;
			 dm_dt[j] = 0;
			 dm_dt[k] = 0;
			 continue;
		}

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

