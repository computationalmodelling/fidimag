#include "nebm_lib.h"
#include "math.h"
#include <stdlib.h>

void cross_product(double *restrict output, double * A, double * B){
    /* Cross product between arrays A and B, assuming they
     * have length = 3
     * Every resulting component is stored in the *output array
     */

    output[0] = A[1] * B[2] - A[2] * B[1];
    output[1] = A[2] * B[0] - A[0] * B[2];
    output[2] = A[0] * B[1] - A[1] * B[0];
}

double dot_product(double * A, double * B, int n){
    /* Dot product between arrays A and B , assuming they
     * have length n */

    double dotp = 0;
    for(int i = 0; i < n; i++){
        dotp += A[i] * B[i];
    }

    return dotp;
}

double compute_norm(double *restrict a, int n, int scale) {
    /* Compute the norm of an array *a

     * This function is necessary to either normalise a vector or compute the
     * distance between two images (see below)
     *
     * ARGUMENTS:

     * n        :: length of a
     *
     * scale    :: If scale is zero, we do not scale the norm

     * Notice that, when computing the DISTANCE between two images, Y_i, Y_i+1
     * for example, we just do the difference (Y_i - Y_(i+1)) and then compute
     * its norm (using this function) but RESCALING the length by the number
     * of components, thus we use scale=n
     *
     */

    double norm = 0;

    for(int i = 0; i < n; i++){
        norm += a[i] * a[i];
    }

    if (scale == 0) norm = sqrt(norm);
    else norm = sqrt(norm) / ((double) scale);

    return norm;
}

void normalise(double *restrict a, int n){

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


void normalise_images_C(double *restrict y, int n_images, int n_dofs_image){

    int i;

    // Index where the components of an image start in the *y array,
    int im_idx;

    for(i = 1; i < n_images - 1; i++){
        im_idx = i * n_dofs_image;
        double * y_i = &y[im_idx];

        normalise(y_i, n_dofs_image);
    }
}


/* ------------------------------------------------------------------------- */


void compute_tangents_C(double *restrict tangents, double *restrict y, double *restrict energies,
                        int n_dofs_image, int n_images
                        ) {

    /* Compute the tangents for every degree of freedom of a NEBM band,
     * following the rules from Henkelman and Jonsson [Henkelman et al.,
     * Journal of Chemical Physics 113, 22 (2000)]
     *
     * We assume we have n_images images in the band, which is represented by
     * the *y array. the *y array has the structure

     *      y = [ theta0_0 phi0_0 theta0_1, phi0_1  ... phi0_(n_dofs_images),
     *            theta1_0 phi1_0 theta1_1, phi1_1  ... phi1_(n_dofs_images),
     *            ...
     *            theta(n_images -1)_0 phi(n_images -1)_0  ... phi(n_images - 1)_(n_dofs_images)]

     * where (theta(i)_j, phi(i)_j) are the spherical coordinates of the j-th
     * spin in the i-th image of the band. Thus we have n_dofs_image spins per
     * image. Images at the extremes, i=0,(n_images-1), are kept fixed, thus we
     * do not compute the tangents for them.
     *
     * The *energies array contains the energy of every image, thus its length
     * is n_images.
     *
     * The tangents have the same length of an image, thus the *tangents array
     * has length (n_images * n_dofs_image), just as the *y array, and we keep
     * the first and last n_dofs_image components as zero (they correspond to
     * the extreme images)
     *
     * To compute the tangents, we use the t+ and t- vectors. Denoting the i-th
     * image by Y_i, these vectors are defined by

     *      t+_i = Y_(i+1) - Y_i          t-_i = Y_i - Y_(i-1)

     * Then, when an image is at a saddle point, we have, using the energies
     * of the neighbouring images,:

     *               __
     *      t_i  =   |   t+_i    if  E_(i+1) > E_i > E_(i-1)
     *              <
     *               |_  t-_i    if  E_(i+1) < E_i < E_(i-1)

     * Otherwise, if we have a maximum or a minimum  (E_(i+1) < E_i > E_(i-1)
     * or E_(i+1) > E_i < E_(i-1) respectively), we use an average of the
     * vectors:
     *
                     __
     *      t_i  =   |   t+_i * dE_max + t-_i * dE_min  if  E_(i+1) > E_(i-1)
     *              <
     *               |_  t+_i * dE_min + t-_i * dE_max  if  E_(i+1) < E_(i-1)

     * where dE_max = max( |E_(i+1) - E_i|, |E_i - E_(i-1)|) and
     *       dE_min = min( |E_(i+1) - E_i|, |E_i - E_(i-1)|)
     *
     */

    int i, j;

    // Index where the components of an image start in the *y array,
    int im_idx;
    // And also the previous and next images:
    int next_im_idx, prev_im_idx;

    double *t_plus;
    double *t_minus;
    t_plus  = malloc(n_dofs_image * sizeof(double));
    t_minus = malloc(n_dofs_image * sizeof(double));

    double deltaE_plus, deltaE_minus;
    double deltaE_MAX, deltaE_MIN;

    for(i = 1; i < n_images - 1; i++){

        im_idx = i * (n_dofs_image);
        next_im_idx = (i + 1) * (n_dofs_image);
        prev_im_idx = (i - 1) * (n_dofs_image);

        // Tangents of the i-th image
        double * t = &tangents[im_idx];

        // Compute the t+ and t- vectors for the i-th image of the band, which
        // is given by the difference of the Y_i image with its neighbours
        for(j = 0; j < n_dofs_image; j++){
            t_plus[j]  = y[next_im_idx + j] - y[im_idx + j];
            t_minus[j] = y[im_idx + j]      - y[prev_im_idx + j];
        }

        // Similarly, compute the energy differences
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

    } // Close loop in images
    free(t_plus);
    free(t_minus);
} // Close main function


/* ------------------------------------------------------------------------- */


void compute_spring_force_C(
        double *restrict spring_force,
        double *restrict y,
        double *restrict tangents,
        double *restrict k,
        int n_images,
        int n_dofs_image,
        double (* compute_distance)(double *, double *, int, int *, int),
        int *restrict material,
        int n_dofs_image_material
        ) {

    /* Compute the spring force for every image of an energy band, which is
     * defined in the *y array. The *y array has (n_images * n_dofs_image)
     * elements. The degrees od freedom (spin components) of the i-th image of
     * the band start at the (i * n_dofs_image) position of the array.
     *
     * For the i-th image, the spring force is defined as:

     *      spring_force_i = k * ( | Y_(i+1) - Y_i| - | Y_i - Y_(i-1)| ) * t_i

     * where t_i is the tangent of the i-th image (see the coresponding
     * function for its definition). Thus, spring_force_i has a length of
     * (n_dofs_image). We do not compute the force for the extremal images,
     * i.e. (i=0, (n_images-1))
     *
     * The norm  | . | between two neighbouring images is calculated as an
     * distance, which depends on the coordinate system chosen.
     *
     * The function: (* compute_distance)(double *, double *, int, int *, int)
     * depends on: (Image_1, 
     *              Image_2, 
     *              number of spins per image,
     *              array indicating dofs (spins) in sites WITH material,
     *              number of dofs in sites with material
     *              ) 
     *
     *  - With material we mean Ms > 0 or mu_s > 0
     */

    int i, j;

    // Index where the components of an image start in the *y array,
    int im_idx;
    // And also the previous and next images:
    int next_im_idx, prev_im_idx;

    double dY_plus_norm, dY_minus_norm;

    for(i = 1; i < n_images - 1; i++){

        im_idx = i * (n_dofs_image);
        next_im_idx = (i + 1) * (n_dofs_image);
        prev_im_idx = (i - 1) * (n_dofs_image);

        double * sf = &spring_force[im_idx];
        double * t = &tangents[im_idx];

        // Compute the distances between the i-th image, Y_i, and its
        // neighbours, the (i+1)-th and (i-1)-th images. The distance between
        // two images is just the norm of their difference scaled by the length
        // of the array (see the compute_norm function)
        dY_plus_norm  = compute_distance(&y[next_im_idx], &y[im_idx],
                                         n_dofs_image,
                                         material, n_dofs_image_material
                                         );
        dY_minus_norm = compute_distance(&y[im_idx], &y[prev_im_idx],
                                         n_dofs_image,
                                         material, n_dofs_image_material
                                         );

        // Now compute the spring force
        for(j = 0; j < n_dofs_image; j++) {
            sf[j] = k[i] * (dY_plus_norm - dY_minus_norm) * t[j];
        }
    }

}


/* ------------------------------------------------------------------------- */


void compute_effective_force_C(double *restrict G,
                               double *restrict tangents,
                               double *restrict gradientE,
                               double *restrict spring_force,
                               int    *restrict climbing_image,
                               int n_images,
                               int n_dofs_image) {

    /* Compute the effective force for every image of an energy band, which is
     * defined in the *y array. The *y array has (n_images * n_dofs_image)
     * elements. The degrees of freedom (spin components) of the i-th image of
     * the band start at the (i * n_dofs_image) position of the array.
     *
     * For the i-th image, the effective force G, is defined as:

     *      G_i = - gradient( E_i ) +  [ gradient(E_i) . t_i ] t_i
     *            + spring_force_i

     * where t_i is the tangent of the i-th image (see the coresponding
     * function for its definition) and e_i is the energy of the i-th image.
     * The gradient is with respect to the degrees of freedom (or magnetisation
     * field). Since, physically, the gradient of the energy is the effective
     * field, we use this fact to compute the gradient, but this is already
     * calculated in the gradientE array.
     *
     * For climbing images we invert the force component and remove the spring
     * force
     *
     */

    int i, j;

    // Index where the components of an image start in the *y array,
    int im_idx;
    double gradE_dot_t;

    for(i = 1; i < n_images - 1; i++){

        im_idx = i * (n_dofs_image);

        double * t = &tangents[im_idx];
        double * gradE = &gradientE[im_idx];
        double * sf = &spring_force[im_idx];

        gradE_dot_t = dot_product(gradE, t, n_dofs_image);

        if(climbing_image[i] == 0) {
            for(j = 0; j < n_dofs_image; j++) {
                G[im_idx + j] = -gradE[j] + gradE_dot_t * t[j] + sf[j];
            }
        }
        else {
            for(j = 0; j < n_dofs_image; j++) {
                G[im_idx + j] = -gradE[j] + 2 * gradE_dot_t * t[j];
            }
        }
    }
}

