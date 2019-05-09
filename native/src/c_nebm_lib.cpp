#include "c_nebm_lib.h"
#include "math.h"
#include <stdlib.h>

void cross_product(double * output, double * A, double * B){
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

double compute_norm(double * a, int n) {
    /* Compute the norm of an array *a
     *
     * ARGUMENTS:

     * n        :: length of a

     */

    double norm = 0;
    for(int i = 0; i < n; i++){ norm += a[i] * a[i]; }
    norm = sqrt(norm);
    return norm;
}

void normalise(double * a, int n){

    /* Normalise the *a array, whose length is n (3 * number of nodes in
     * cartesian, and 2 * number of nodes in spherical) To do this we compute
     * the length of *a :
     *
     *      SQRT[ a[0] ** 2 + a[1] ** 2 + ... ]
     *
     *  and divide every *a entry by that length
     */

    double length;

    length = compute_norm(a, n);

    if (length > 0){
        length = 1.0 / length;
    }

    for(int i = 0; i < n; i++){
        a[i] *= length;
    }
}


void normalise_images_C(double * y, int n_images, int n_dofs_image){

    int i;

    // Index where the components of an image start in the *y array,
    int im_idx;

    for(i = 1; i < n_images - 1; i++){
        im_idx = i * n_dofs_image;
        double * y_i = &y[im_idx];

        normalise(y_i, n_dofs_image);
    }
}

void normalise_spins_C(double * y, int n_images, int n_dofs_image){

    int i, j;

    // Index where the components of an image start in the *y array,
    int im_idx;

    for(i = 1; i < n_images - 1; i++){
        im_idx = i * n_dofs_image;

        // Normalise per spin
        for(j = 0; j < n_dofs_image; j+=3){
            double * y_i = &y[im_idx + j];
            normalise(y_i, 3);
        }
    }
}


/* ------------------------------------------------------------------------- */


void compute_tangents_C(double * tangents, double * y, double * energies,
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
    t_plus  = new double[n_dofs_image];
    t_minus = new double[n_dofs_image];

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

        else if((deltaE_plus < 0 && deltaE_minus > 0) ||      // A maximum
                (deltaE_plus > 0 && deltaE_minus < 0)) {      // A minimum

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

        // When energy has a constant slope
        else {
            for(j = 0; j < n_dofs_image; j++) t[j] = t_plus[j] + t_minus[j];
        }

        /* ----------------------------------------------------------------- */

    } // Close loop in images
    delete[] t_plus;
    delete[] t_minus;
} // Close main function


/* ------------------------------------------------------------------------- */

void compute_spring_force_C(
        double * spring_force,
        double * y,
        double * tangents,
        double * k,
        int n_images,
        int n_dofs_image,
        double * distances
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
     * The norm  | . | between two neighbouring images is calculated as a
     * distance, which depends on the chosen coordinate system. Distances
     * are stored in the *distances array which has the differences
     *      [ |Y0 - Y1|, |Y1 - Y2|, |Y2 - Y3|, ...]
     *
     */

    int i, j;

    // Index where the components of an image start in the *y array,
    int im_idx;

    double dY_plus_norm, dY_minus_norm;

    for(i = 1; i < n_images - 1; i++){

        im_idx = i * (n_dofs_image);

        double * sf = &spring_force[im_idx];
        double * t = &tangents[im_idx];

        // Get the distances between the i-th image, Y_i, and its
        // neighbours, the (i+1)-th and (i-1)-th images.
        dY_plus_norm  = distances[i];
        dY_minus_norm = distances[i - 1];

        // Now compute the spring force
        for(j = 0; j < n_dofs_image; j++) {
            sf[j] = k[i] * (dY_plus_norm - dY_minus_norm) * t[j];
        }
    }

}


/* ------------------------------------------------------------------------- */


void compute_effective_force_C(double * G,
                               double * tangents,
                               double * gradientE,
                               double * spring_force,
                               int    * climbing_image,
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
        // falling image with spring constant -> 0
        else if(climbing_image[i] == -1) {
            for(j = 0; j < n_dofs_image; j++) {
                G[im_idx + j] = -gradE[j];
            }
        }
        // climbing image
        else {
            for(j = 0; j < n_dofs_image; j++) {
                G[im_idx + j] = -gradE[j] + 2 * gradE_dot_t * t[j];
            }
        }
    }
}

// ----------------------------------------------------------------------------

// Compute the distances between consecutive images in the band, using the compute_distance
// function that depends on the chosen coordinate system to represent the
// degrees of freedom. This function updates the *distances and *path_distances
// arrays. The *distances array saves the distances as:
//      [ |Y1 - Y0|, |Y2 - Y1|, ... ]
// Path distances are the total distances relative to the 0-th image 
void compute_image_distances(double * distances,
                             double * path_distances,
                             double * y,
                             int n_images,
                             int n_dofs_image,
                             double (* compute_distance)(double *,
                                                         double *,
                                                         int,
                                                         int *,
                                                         int),
                             int *  material,
                             int n_dofs_image_material
                             ) {

    int i, im_idx, next_im_idx;
    path_distances[0] = 0.0;
    for(i = 0; i < n_images - 1; i++){

        im_idx = i * (n_dofs_image);
        next_im_idx = (i + 1) * (n_dofs_image);

        distances[i] = compute_distance(&y[next_im_idx], &y[im_idx],
                                        n_dofs_image,
                                        material,
                                        n_dofs_image_material
                                        );
        path_distances[i + 1] = path_distances[i] + distances[i];
    }
}

// ----------------------------------------------------------------------------

void project_vector_C(double * vector, double * y_i,
                      int n_dofs_image
                      ){

    /* Project a vector into the space of the y_i image in an energy band. We
     * assume that *image and *y_i have the same length. The vector and the
     * image are in Cartesian coordinates. Thus, remembering that an image y_i
     * has the components:

     *      [ spin0_x spin0_y spin0_z spin1_x ...]

     * and assuming that *vector has the same order, the projection is made
     * through the components of every 3 component vector in the arrays. This
     * means, if V_i is the i-th vector of the *vector array and m_i is the
     * i-th spin in the y_i array, we have to compute

     *                  ->      ->      ->    ->   ->
     *      Projection( V_i ) = V_i - ( V_i o m_i) m_i
     *                               (dot product)

     * for every V_i in vector. The components of the V_i vector start at the
     * (3 * i) position of the *vector array, and the i-th spin components
     * start at the (3 * i) position of the y_i array
     *
     * - Remember that an image has n_dofs_image elements, thus the number of
     *   spins is n_dofs_image / 3
     *
     */

    int j;
    double v_dot_m_i = 0;

    for(j = 0; j < n_dofs_image; j++) {
        // Every 3 components of the *vector and *y_i array, compute the
        // dot product of the i-th vector (i.e. using
        //      ( vector[j], vector[j+1], vector[j+2] ) , j=0,3,6, ... )
        // and then compute the projection for the 3 components of the vector
        if (j % 3 == 0) v_dot_m_i = dot_product(&vector[j], &y_i[j], 3);
        vector[j] = vector[j] - v_dot_m_i * y_i[j];
    }
}


void project_images_C(double * vector, double * y,
                      int n_images, int n_dofs_image
                      ){

    /* Compute the projections of the vectors in the *vector array into the
     * space formed by the spin vectors in the array *y
     *
     * We assume that the *vector array is made of n_images images, i.e.

     *      vector = [ vector(0)0_x vector(0)0_y vector(0)0_z vector(0)1_x ... ] IMAGE 0
                       vector(1)0_x vector(1)0_y vector(1)0_z vector(1)1_x ... ] IMAGE 1
                       ...                                                        ...
     *                ]

     * where vector(i)j_x is the x component of the j-th vector of the i-th
     * image (similar for y_i, only that vectors are spins). Thus, every image
     * starts at the i * n_dofs_image position of the *vector and *y_i array
     *
     * IMPORTANT: projections are NOT calculated for the extrema images,
     *            i.e. for IMAGE_0 and IMAGE_N
     *
     * - Notice that we could have just computed the projections using the
     *   whole vector array, without separating it into images, but the
     *   approach taken here is clearer if we use the project_vector_C
     *   function. In addition, we can parallelise the calculation in the
     *   future, computing the projection of images in different threads.
     *
     * - An image is simply a copy of the magnetic system, i.e. every image has
     *   n_dofs_image elements
     */

    int i;

    // Index where the components of an image start in the *y array,
    int im_idx;

    for(i = 1; i < n_images - 1; i++){

        im_idx = i * (n_dofs_image);

        double * v = &vector[im_idx];
        double * y_i = &y[im_idx];

        project_vector_C(v, y_i, n_dofs_image);

    }
}

// ----------------------------------------------------------------------------

double compute_distance_cartesian(double * A, double * B, int n_dofs_image,
                                  int * material, int n_dofs_image_material
                                  ) {

    /* Compute the distance between two images, A and B, discarding the sites
     * without material
     *
     * We still consider pinned sites because they contribute to the scale
     * factor for the distance
     */

    double A_minus_B[n_dofs_image_material];
    double distance;
    int j = 0;

    for(int i = 0; i < n_dofs_image; i++) {
        if (material[i] > 0) {
            A_minus_B[j]  = A[i] - B[i];
            j += 1;
        }
    }

    // CHECK division by number of dofs
    distance = compute_norm(A_minus_B, n_dofs_image_material) / n_dofs_image_material;

    return distance;
}

// ----------------------------------------------------------------------------

/* Landau-Lifshitz like equation (no precession term) for the evolution of the
 * NEBM in Cartesian coordinates. To preserve the magnitude of the spins we add
 * a correction term. Calling G_eff to the effective force, every image Y of a
 * NEBM band evolves according to the equation:

 *      dY                                     2
 *     ---- = - Y x ( Y x G_eff ) + c ( 1 - Y  ) Y
 *      dt

 *  with
                    _________
 *                  /        2
 *      c = 6 *    / (  dY  )
 *                /  ( ---- )
 *              \/   (  dt  )
 *
 */

void compute_dYdt(double * Y, double * G, double * dYdt,
                  int * pins,
                  int n_dofs_image
                  ) {

    int n_spins = n_dofs_image / 3;
    for(int i = 0; i < n_spins; i++){
       	int j = 3 * i;

        if (pins[i] > 0){
            dYdt[j] = 0;
		    dYdt[j + 1] = 0;
		    dYdt[j + 2] = 0;
		    continue;
		}

        double Y_dot_Y = Y[j] * Y[j] + Y[j + 1] * Y[j + 1] + Y[j + 2] * Y[j + 2];
       	double Y_dot_G = Y[j] * G[j] + Y[j + 1] * G[j + 1] + Y[j + 2] * G[j + 2];
        // (Y * Y) G - (Y * G) Y = - Y x (Y x G)
       	dYdt[j] = Y_dot_Y * G[j] - Y_dot_G * Y[j];
       	dYdt[j + 1] = Y_dot_Y * G[j + 1] - Y_dot_G * Y[j + 1];
       	dYdt[j + 2] = Y_dot_Y * G[j + 2] - Y_dot_G * Y[j + 2];

        // Correction factor to rescale the spin length at every iteration step
       	double c = 6 * sqrt(dYdt[j]     * dYdt[j]     +
                            dYdt[j + 1] * dYdt[j + 1] +
                            dYdt[j + 2] * dYdt[j + 2]);

       	dYdt[j]     += c * (1 - Y_dot_Y) * Y[j];
        dYdt[j + 1] += c * (1 - Y_dot_Y) * Y[j + 1];
       	dYdt[j + 2] += c * (1 - Y_dot_Y) * Y[j + 2];
    }
}

void compute_dYdt_C(double * y, double * G, double * dYdt, int * pins,
                    int n_images, int n_dofs_image) {
    #pragma omp parallel for schedule(static)
	for(int i = 1; i < n_images - 1; i++){
        //printf("");
        int j = i * n_dofs_image;
        compute_dYdt(&y[j], &G[j], &dYdt[j], &pins[0], n_dofs_image);
    }
    return;
}

// ----------------------------------------------------------------------------
// dYdt functions without the magic correction factor
// Necessary to implement more standard integrators, such as Verlet

void compute_dYdt_nc(double * Y, double * G, double * dYdt,
                     int * pins,
                     int n_dofs_image
                     ) {

    int n_spins = n_dofs_image / 3;
    for(int i = 0; i < n_spins; i++){
       	int j = 3 * i;

        if (pins[i] > 0){
            dYdt[j] = 0;
		    dYdt[j + 1] = 0;
		    dYdt[j + 2] = 0;
		    continue;
		}

        double Y_dot_Y = Y[j] * Y[j] + Y[j + 1] * Y[j + 1] + Y[j + 2] * Y[j + 2];
       	double Y_dot_G = Y[j] * G[j] + Y[j + 1] * G[j + 1] + Y[j + 2] * G[j + 2];
        // (Y * Y) G - (Y * G) Y = - Y x (Y x G)
       	dYdt[j] = Y_dot_Y * G[j] - Y_dot_G * Y[j];
       	dYdt[j + 1] = Y_dot_Y * G[j + 1] - Y_dot_G * Y[j + 1];
       	dYdt[j + 2] = Y_dot_Y * G[j + 2] - Y_dot_G * Y[j + 2];
    }
}

void compute_dYdt_nc_C(double * y, double * G, double * dYdt,
                       int * pins,
                       int n_images, int n_dofs_image) {
    #pragma omp parallel for schedule(static)
	for(int i = 1; i < n_images - 1; i++){
        //printf("");
        int j = i * n_dofs_image;
        compute_dYdt_nc(&y[j], &G[j], &dYdt[j], &pins[0], n_dofs_image);
    }
    return;
}
