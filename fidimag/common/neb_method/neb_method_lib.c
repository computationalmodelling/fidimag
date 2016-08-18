#include "neb_method_lib.h"
#include "math.h"

inline double compute_norm(double *a, int n, int scale) {
    /* a        :: array
     * n        :: length of a
     * scale    :: If scale is zero, we do not scale the norm
     */

    double norm = 0;
    double scale_db = (double) scale;

    for(int i = 0; i < n; i++){

        if (a[i] > WIDE_PI){
            a[i] = 2 * WIDE_PI - a[i];
        } else if(a[i] < -WIDE_PI){
            a[i] += 2 * WIDE_PI;
        }

        norm += a[i] * a[i];
    }

    if (scale == 0) norm = sqrt(norm);
    else norm = sqrt(norm) / scale_db;
    
    return norm;
}

inline void normalise(double *a, int n){

    /* Normalise the *a array, whose length is n (3 * number of nodes in
     * cartesian, and 2 * number of nodes in spherical) To do this we compute
     * the length of *a :
     *
     *      SQRT[ a[0] ** 2 + a[1] ** 2 + ... ]
     *
     *  and divide every *a entry by that length
     */

    double length;

    length = compute_norm(a, n, 0);

    if (length > 0){
        length = 1.0 / length;
    }

    for(int i = 0; i < n; i++){
        a[i] *= length;
    }
}

void compute_tangents_C(double *tangents, double *y, double *energies,
                        int n_dofs_image, int n_images
                        ) {

    int i, j;

    // Index where the components of an image start in the *y array,
    int im_idx;
    // And also the previous and next images:
    int next_im_idx, prev_im_idx;

    double t_plus[n_dofs_image];
    double t_minus[n_dofs_image];

    double deltaE_plus, deltaE_minus;
    double deltaE_MAX, deltaE_MIN;

    for(i = 1; i < n_images -1; i++){

        im_idx = i * (n_dofs_image);
        next_im_idx = (i + 1) * (n_dofs_image);
        prev_im_idx = (i - 1) * (n_dofs_image);

        // Tangents of the i-th image
        double *t = &tangents[im_idx];

        for(j = 0; j < n_dofs_image; j++){
            t_plus[j]  = y[next_im_idx + j] - y[im_idx + j];
            t_minus[j] = y[im_idx + j]      - y[prev_im_idx + j];
        }

        deltaE_plus  = energies[i + 1] - energies[i];
        deltaE_minus = energies[i]     - energies[i - 1];

        /* Now we follow Henkelman and Jonsson rules [Henkelman et al., Journal
         * of Chemical Physics 113, 22 (2000)] for the tangent directions
         * (remember that there is a tangent for every spin (degree of
         * freedom))
         *
         * The first two cases are straightforward: If the energy has a
         * positive (negative) slope, just make the difference between the spin
         * directions with respect to the right (left) image.
         *
         * The other case is when the image is an energy maximum, or minimum
         */
        
        /* ----------------------------------------------------------------- */

        if(deltaE_plus > 0 && deltaE_minus > 0) {
            for(j = 0; j < n_dofs_image; j++) t[j] = t_plus[j];
        }

        /* ----------------------------------------------------------------- */

        else if(deltaE_plus < 0 && deltaE_minus < 0) {
            for(j = 0; j < n_dofs_image; j++) t[j] = t_minus[j];
        }

        /* ----------------------------------------------------------------- */

        else {
            
            /* According to the energy of the neighbouring images, the tangent
             * of the i-th image will be a combination of the differences wrt
             * to the left and right images components weighted according to
             * the neighbours energies
             */
            deltaE_MAX = fmax(fabs(deltaE_plus), fabs(deltaE_minus));
            deltaE_MIN = fmin(fabs(deltaE_plus), fabs(deltaE_minus));

            if (energies[i + 1] > energies[i - 1]) {
                for(j = 0; j < n_dofs_image; j++) {
                    t[j] = t_plus[j] * deltaE_MAX + t_minus[j] * deltaE_MIN;
                }
            }

            else {
                for(j = 0; j < n_dofs_image; j++) {
                    t[j] = t_plus[j] * deltaE_MIN + t_minus[j] * deltaE_MAX;
                }
            }

        }

        /* ----------------------------------------------------------------- */

        normalise(t, n_dofs_image);

    } // Close loop in images

} // Close main function


/* ------------------------------------------------------------------------- */

// void project_tangents(double *tangents, double *y){
// 
// }

/* ------------------------------------------------------------------------- */

void compute_spring_force_C(double *spring_force,
                            double *y,
                            double *tangents,
                            double k,
                            int n_images,
                            int n_dofs_image
                            ) {

    int i, j;

    // Index where the components of an image start in the *y array,
    int im_idx;
    // And also the previous and next images:
    int next_im_idx, prev_im_idx;

    double dY_plus[n_dofs_image];
    double dY_minus[n_dofs_image];

    double dY_plus_norm, dY_minus_norm;

    for(i = 1; i < n_images -1; i++){

        im_idx = i * (n_dofs_image);
        next_im_idx = (i + 1) * (n_dofs_image);
        prev_im_idx = (i - 1) * (n_dofs_image);

        double * sf = &spring_force[im_idx];
        double * t = &tangents[im_idx];

        // Compute the distances between the i-th image, Y_i, and its
        // neighbours, the (i+1)-th and (i-1)-th images. The distance between
        // two images is just the norm of their difference scaled by the length
        // of the array
        for(j = 0; j < n_dofs_image; j++) {
            dY_plus[j]  = y[next_im_idx + j] - y[im_idx + j];
            dY_minus[j] = y[im_idx + j]      - y[prev_im_idx + j];
        }
        dY_plus_norm  = compute_norm(dY_plus,  n_dofs_image, n_dofs_image);
        dY_minus_norm = compute_norm(dY_minus, n_dofs_image, n_dofs_image);

        // Now compute the spring force
        for(j = 0; j < n_dofs_image; j++) {
            sf[j] = k * (dY_plus_norm - dY_minus_norm) * t[j];
        }
    }

}
