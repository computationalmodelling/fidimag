#include "clib.h"

//compute (\vec{j}_s \cdot \nabla \vec{m})
void compute_stt_field_c(double *spin, double *field, double jx, double jy, double jz,
                       double dx, double dy, double dz, int nx, int ny, int nz) {
    
	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int id, id1, id2;
	
    
    for(i=0;i<3*n1;i++){
        field[i]=0;
    }
    
    if (jx!=0 && nx>1){
        
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                for (k = 0; k < nz; k++) {
                    id = nyz * i + nz * j + k;
                    
                    id1 = id - nyz;
                    id2 = id + nyz;
                    
                    if (i==0){
                        id1 = id;
                    }else if (i == nx-1){
                        id2 = id;
                    }
                    
		    field[id] += jx * (spin[id2]-spin[id1])/(2*dx);
                    field[id+n1] += jx * (spin[id2+n1]-spin[id1+n1])/(2*dx);
                    field[id+n2] += jx * (spin[id2+n2]-spin[id1+n2])/(2*dx);
                }
            }
        }
    }
    
    
    if (jy!=0 && ny>1){
        
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                for (k = 0; k < nz; k++) {
                    id = nyz * i + nz * j + k;
                    
                    id1 = id - nz;
                    id2 = id + nz;
                    
                    if (j==0){
                        id1 = id;
                    }else if (j == ny-1){
                        id2 = id;
                    }
                    
		    field[id] += jy * (spin[id2]-spin[id1])/(2*dy);
                    field[id+n1] += jy * (spin[id2+n1]-spin[id1+n1])/(2*dy);
                    field[id+n2] += jy * (spin[id2+n2]-spin[id1+n2])/(2*dy);
                }
            }
        }
    }
    
    if (jz!=0 && nz>1){
        
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                for (k = 0; k < nz; k++) {
                    id = nyz * i + nz * j + k;
                    
                    id1 = id - 1;
                    id2 = id + 1;
                    
                    if (k==0){
                        id1 = id;
                    }else if (k == nz-1){
                        id2 = id;
                    }
                    
	            field[id] += jz * (spin[id2]-spin[id1])/(2*dz);
                    field[id+n1] += jz * (spin[id2+n1]-spin[id1+n1])/(2*dz);
                    field[id+n2] += jz * (spin[id2+n2]-spin[id1+n2])/(2*dz);
                }
            }
        }
    }
    
}



void llg_stt_rhs(double *dm_dt, double *m, double *h, double *h_stt, double *alpha,
                   double beta, double u0, double gamma, int nxyz) {
    
	int i, j, k;
    
	double mth0, mth1, mth2;
    double coeff, coeff_stt, mh, mht;
    double hpi,hpj,hpk;
    
	for (i = 0; i < nxyz; i++) {
		j = i + nxyz;
		k = j + nxyz;
        
		coeff = -gamma/(1+alpha[i]*alpha[i]);
        
        mh = m[i]*h[i] + m[j]*h[j] + m[k]*h[k];
        
        //we use the same trick as llg does that drops the mm
        hpi = h[i] - mh*m[i];
        hpj = h[j] - mh*m[j];
        hpk = h[k] - mh*m[k];
        
        mth0 = (m[j] * hpk - m[k] * hpj);
		mth1 = (m[k] * hpi - m[i] * hpk);
		mth2 = (m[i] * hpj - m[j] * hpi);
        
        dm_dt[i] = coeff*(mth0 - hpi * alpha[i]);
        dm_dt[j] = coeff*(mth1 - hpj * alpha[i]);
        dm_dt[k] = coeff*(mth2 - hpk * alpha[i]);
        
        //the above part is standard LLG equation.
        
        coeff_stt = u0/(1+alpha[i]*alpha[i]);
        
        mht = m[i]*h_stt[i] + m[j]*h_stt[j] + m[k]*h_stt[k];
        
        hpi = h_stt[i] - mht*m[i];
        hpj = h_stt[j] - mht*m[j];
        hpk = h_stt[k] - mht*m[k];
        
        mth0 = (m[j] * hpk - m[k] * hpj);
        mth1 = (m[k] * hpi - m[i] * hpk);
        mth2 = (m[i] * hpj - m[j] * hpi);
        
        dm_dt[i] += coeff_stt*((1 + alpha[i]*beta)*hpi - (beta - alpha[i])*mth0);
        dm_dt[j] += coeff_stt*((1 + alpha[i]*beta)*hpj - (beta - alpha[i])*mth1);
        dm_dt[k] += coeff_stt*((1 + alpha[i]*beta)*hpk - (beta - alpha[i])*mth2);
        
        
	}
    
}













