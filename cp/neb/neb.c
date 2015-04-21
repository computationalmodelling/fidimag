#include "neb.h"


//normalise the given array
inline void normalise(double *res, int n){

		double length=0;

		for(int i=0;i<n;i++){
			length += res[i]*res[i];
		}

		if (length>0){
			length = 1.0/sqrt(length);
		}

		for(int i=0;i<n;i++){
			res[i]*=length;
		}

}

//compute b-a and store it to res with normalised length 1
inline void difference(double *res, double *a, double *b, int n){

		for(int i=0;i<n;i++){
			res[i] = b[i]-a[i];
		}

		normalise(res, n);

}


void compute_tangents_c(double *ys, double *energy, double *tangents, int total_img_num, int nodes) {
	//length of ys is total_img_num * nodes
	//length of energy is total_img_num-2
	//length of tangents is (total_img_num-2) * nodes
	//nodes = 3*nxyz

	double t1[nodes];
	double t2[nodes];

	for(int i=1; i<total_img_num-1; i++){
		int ja = (i-1)*nodes;
		int jb = i*nodes;
		int jc = (i+1)*nodes;
	   	double *ya = &ys[ja];
	   	double *yb = &ys[jb];
	   	double *yc = &ys[jc];
	            	
	    double *t = &tangents[ja];

	    double e1 = energy[i-1]-energy[i];
	    double e2 = energy[i]-energy[i+1];

	    if (e1<0 && e2<0){
	            difference(t, yb, yc, nodes);
	    }else if(e1>0&&e2>0){
	            difference(t, ya, yb, nodes);
	    }else{
	            		difference(t1, ya, yb, nodes);
	            		difference(t2, yb, yc, nodes);

	            		double max_e, min_e;

	            		if (fabs(e1)>fabs(e2)){
	            			max_e = fabs(e1);
	            			min_e = fabs(e2);
	            		}else{
	            			max_e = fabs(e2);
	            			min_e = fabs(e1);
	            		}

	            		if (energy[i+1]>energy[i-1]){
	            			for(int i=0;i<nodes;i++){
	            				t[i] = min_e*t1[i] + max_e*t2[i];
	            			}
	            		}else{
	            			for(int i=0;i<nodes;i++){
	            				t[i] = max_e*t1[i] + min_e*t2[i];
	            			}
	            		}

	            		normalise(t, nodes);
	            	}


	            }

		}


inline void compute_dmdt(double *m, double *h, double *dm_dt, int *pins, int nodes){

    int n = nodes/3;

    for(int i=0;i<n;i++){
       	int j=i+n;
        int k=j+n;

        if (pins[i]>0){
			 dm_dt[i] = 0;
			 dm_dt[j] = 0;
			 dm_dt[k] = 0;
			 continue;
		}

        double mm = m[i]*m[i] + m[j]*m[j] + m[k]*m[k];
       	double mh = m[i]*h[i] + m[j]*h[j] + m[k]*h[k];
        //mm.h-mh.m=-mx(mxh)
       	dm_dt[i] = mm*h[i] - mh*m[i];
       	dm_dt[j] = mm*h[j] - mh*m[j];
       	dm_dt[k] = mm*h[k] - mh*m[k];

       	double c = 6*sqrt(dm_dt[i]*dm_dt[i]+dm_dt[j]*dm_dt[j]+dm_dt[k]*dm_dt[k]);

       	dm_dt[i] += c*(1-mm)*m[i];
        dm_dt[j] += c*(1-mm)*m[j];
       	dm_dt[k] += c*(1-mm)*m[k];
    }


}

void compute_dm_dt_c(double *ys, double *heff, double *dm_dt, int *pins, int image_num, int nodes) {

	for(int i=1; i<image_num-1; i++){

		int j = i*nodes;
		
		compute_dmdt(&ys[j], &heff[j], &dm_dt[j], pins, nodes);

        }

    return;

}
