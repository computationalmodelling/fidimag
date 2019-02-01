#include<stdio.h>
#include "nebm_cartesian_lib.h"
#include "nebm_lib.h"
#include "math.h"


void project_vector_C(double *restrict vector, double *restrict y_i,
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


void project_images_C(double *restrict vector, double *restrict y,
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


double compute_distance_cartesian(double *restrict A, double *restrict B, int n_dofs_image,
                                  int *restrict material, int n_dofs_image_material
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

    distance = compute_norm(A_minus_B, n_dofs_image_material,
                            n_dofs_image_material);

    return distance;
}

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

void compute_dYdt(double *restrict Y, double *restrict G, double *restrict dYdt,
                  int *restrict pins,
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

void compute_dYdt_C(double *restrict y, double *restrict G, double *restrict dYdt, int *restrict pins, 
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

void compute_dYdt_nc(double *restrict Y, double *restrict G, double *restrict dYdt,
                     int *restrict pins,
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

void compute_dYdt_nc_C(double *restrict y, double *restrict G, double *restrict dYdt,
                       int *restrict pins, 
                       int n_images, int n_dofs_image) {
    #pragma omp parallel for schedule(static)
	for(int i = 1; i < n_images - 1; i++){
        //printf("");
        int j = i * n_dofs_image;
        compute_dYdt_nc(&y[j], &G[j], &dYdt[j], &pins[0], n_dofs_image);
    }
    return;
}
