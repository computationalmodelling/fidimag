#include "neb.h"
#define PI 3.14159265358979323846

//normalise the given array
inline void normalise(double *res, int n){
    /* Normalise the *res array, whose length is n
     * (3 * number of nodes in cartesian, and
     *  2 * number of nodes in spherical)
     * To do this we compute the length of res:
     *
     *      SQRT[ res[0] ** 2 + res[1] ** 2 + ... ]
     *
     *  and divide every res entry by that length
     */
    double length=0;

    for(int i = 0; i < n; i++){
        
        if (res[i] > PI){
            res[i] = 2 * PI - res[i];
        } else if(res[i] < -PI){
            res[i] += 2 * PI;
        }
        length += res[i] * res[i];
    }

    if (length > 0){
        length = 1.0 / sqrt(length);
    }

    for(int i = 0; i < n; i++){
        res[i] *= length;
    }
}

// compute b-a and store it to res with normalised length 1
inline void difference(double *res, double *a, double *b, int n){

    for(int i = 0; i < n; i++){
        res[i] = b[i] - a[i];
    }
    normalise(res, n);
}


void compute_tangents_c(double *ys, double *energy, double *tangents,
                        int total_img_num, int nodes) {
	/* We have the variables:
     *
     * length of ys is total_img_num * nodes
	 * length of energy is total_img_num-2  ??
	 * length of tangents is (total_img_num-2) * nodes
	 * nodes = 3*nxyz (Cartesian) or 2*nxyz (spherical)
     *
     * The order of the spin directions does not matter, as long as it
     * is the same for every image, i.e. we can use
     *   [mx1, my1, mz1, mx2 ...]  OR
     *   [mx1, mx2, ... my1, my2... , mz1, mz2 ,... ]
     */


	double t1[nodes];
	double t2[nodes];

    /* We will set the tangents for every image. A tangent array has
     * length equal to 3 * nxyz (x,y,z components for every
     * degree of freedom) */
	for(int i=1; i<total_img_num-1; i++){
		/* Indexes for the starting memory location
         * of nearest neighbouring images
         *      ja --> image - 1
         *      jb --> image
         *      jc --> image + 1
         * */
        int ja = (i - 1) * nodes;
		int jb = i * nodes;
		int jc = (i + 1) * nodes;

        /* Starting memory locations of the NN images spin components */
	   	double *ya = &ys[ja];  // y(i - 1)
	   	double *yb = &ys[jb];  // y(i)
	   	double *yc = &ys[jc];  // y(i + 1)

        /* Starting location for the tangent components
         * of the (i - 1) image */
	    double *t = &tangents[ja];

        /* e1 --> energy difference with the neighbouring image to the left
         * e2 --> energy difference with the neighbouring image to the right
         */
	    double e1 = energy[i - 1] - energy[i];
	    double e2 = energy[i] - energy[i + 1];

        /* Now we follow Henkelman and Jonsson rules [Henkelman et al., Journal
         * of Chemical Physics 113, 22 (2000)] for the tangent directions
         * (remember that there is a tangent for every node (degree of freedom,
         * spin) and every node has 3 components)
         *
         * The first two cases are straightforward: If the energy has a
         * positive (negative) slope, just make the difference between the
         * spin directions with respect to the right (left) image.
         *
         * The other case is when the image is an energy saddle point,
         * maximum, or minimum
         */

	    if (e1 < 0 && e2 < 0) {
            /* y(i+ 1)  - y(i) */
	        difference(t, yb, yc, nodes);
	    } else if (e1 > 0 && e2 > 0) {
	        /* y(i)  - y(i - 1)*/
            difference(t, ya, yb, nodes);
	    } else {
            /* t1 --> spin components difference with the left image
             *        i.e. y(i) - y(i - 1)
             * t2 --> spin components difference with the right image
             *        i.e. y(i + 1) - y(i)
             *
             * The difference function saves the results in the t1 and t2
             * arrays
             */
            difference(t1, ya, yb, nodes);
            difference(t2, yb, yc, nodes);

            double max_e, min_e;

            /* Compare energy differences (wrt left or right image)
             * magnitudes */
            if (fabs(e1) > fabs(e2)){
                max_e = fabs(e1);
                min_e = fabs(e2);
            } else {
                max_e = fabs(e2);
                min_e = fabs(e1);
            }

            /* According to the energy of the neighbouring images, the tangent
             * of the i-th image will be a combination of the differences wrt
             * to the left and right images (see above for t1, t2), with their
             * components weighted according to the neighbours energies
             */
            if (energy[i + 1] > energy[i - 1]){
                for(int i=0; i < nodes; i++){
                    t[i] = min_e * t1[i] + max_e * t2[i];
                }
            } else {
                for(int i = 0; i < nodes; i++){
                    t[i] = max_e * t1[i] + min_e * t2[i];
                }
            }
            /* Normalise as in the differences function */
            normalise(t, nodes);

        } // close else
    } // close for loop for images
} // close main function


inline void compute_dmdt(double *m, double *h, double *dm_dt,
                         int *pins, int nodes){

    for(int i=0; i < nodes; i++){
       	int j = 3 * i;

        if (pins[i]>0){
             dm_dt[j] = 0;
			 dm_dt[j + 1] = 0;
			 dm_dt[j + 2] = 0;
			 continue;
		}

        double mm = m[j] * m[j] + m[j + 1] * m[j + 1] + m[j + 2] * m[j + 2];
       	double mh = m[j] * h[j] + m[j + 1] * h[j + 1] + m[j + 2] * h[j + 2];
        //mm.h-mh.m=-mx(mxh)
       	dm_dt[j] = mm * h[j] - mh * m[j];
       	dm_dt[j + 1] = mm * h[j + 1] - mh * m[j + 1];
       	dm_dt[j + 2] = mm * h[j + 2] - mh * m[j + 2];

       	double c = 6 * sqrt(dm_dt[j]     * dm_dt[j]     +
                            dm_dt[j + 1] * dm_dt[j + 1] +
                            dm_dt[j + 2] * dm_dt[j + 2]);

       	dm_dt[j]     += c * (1 - mm) * m[j];
        dm_dt[j + 1] += c * (1 - mm) * m[j + 1];
       	dm_dt[j + 2] += c * (1 - mm) * m[j + 2];
    }

}

void compute_dm_dt_c(double *ys, double *heff, double *dm_dt,
                     int *pins, int image_num, int nodes) {

	for(int i = 1; i < image_num-1; i++){

		int j = 3*i*nodes;

		compute_dmdt(&ys[j], &heff[j], &dm_dt[j], pins, nodes);

        }

    return;

}
