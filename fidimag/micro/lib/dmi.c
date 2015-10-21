#include "micro_clib.h"
#include <stdlib.h>

/* Functions to Compute the micromagnetic Dzyaloshinskii Moriya interaction
*  field and energy using the
* matrix of neighbouring spins and a second order approximation
* for the derivative
*
* Common inputs:
*
* Ms_inv     :: Array with the (1 / Ms) values for every mesh node.
*               The values are zero for points with Ms = 0 (no material)
*
* D          :: Array with the DMI constant values
*
* dx, dy, dz :: Mesh spacings in the corresponding directions
*
* nxyz       :: Number of mesh nodes
*
* ngbs       :: The array of neighbouring spins, which has (6 * nxyz)
*               entries. Specifically, it contains the indexes of
*               the neighbours of every mesh node, in the following order:
*                      -x, +x, -y, +y, -z, +z
*
*               Thus, the array is like:
*              | 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, 1-x, 1+x, 1-y, ...  |
*                i=0                           i=1                ...
*
*              where  0-y  is the index of the neighbour of the 0th spin,
*              in the -y direction, for example. The index value for a
*              neighbour where Ms = 0, is evaluated as -1. The array
*              automatically gives periodic boundaries.
*
*              A basic example is a 3 x 3 two dimensional mesh with PBCs
*              in the X and Y direction:
*
*                                  +-----------+
*                                  | 6 | 7 | 8 |
*                                  +-----------+
*                                  | 3 | 4 | 5 |
*                                  +-----------+
*                                  | 0 | 1 | 2 |
*                                  +-----------+
*
*              so, the first 6 entries (neighbours of the 0th mesh node)
*              of the array would be: [ 2  1  6  3 -1 -1 ... ]
*              (-1 since there is no material in +-z, and a '2' first,
*              since it is the left neighbour which is the PBC in x, etc..)
*
*  For the DMI computation, the field  will depend on the interaction
*  structure. Common types are Interfacial and Bulk DMI. For the
*  Bulk DMI we use the interaction found in materials with
*  point group C_nv (although we can easily generalise the code for
*  other symmetry classes)
*
*  The calculation will be based on the "direction" of the DMI vectors,
*  This will be explained inside the functions
*
*  The *m array contains the spins as:
*      [mx0, my0, mz0, mx1, my1, mz1, mx2, ...]
*  so if we want the starting position of the magnetisation for the
*  i-th spin, we only have to do (3 * i) for mx, (3 * i + 1) for my
*  and (3 * i + 2) for mz
*
*
*  IMPORTANT: The DMI field usually has the structure:
*                   2 * D / (mu0 Ms ) * (Some first  derivative of M)
*       When discretising the derivative, it carries a "2" in the
*       denominator which we "cancel" with the "2" in the prefactor,
*       hence we do not put it explicitly in the calculations
*
*       So, when computing the energy: (-1/2) * mu * Ms * H_dmi
*       we only put the 0.5 factor and don't worry about the "2"s in the
*       field
*
*/

void dmi_field_bulk(double *m, double *field, double *energy, double *Ms_inv,
                    double *D, double dx, double dy, double dz,
                    int nxyz, int *ngbs) {

    /* In the atomic model, the effective field in the i-th spin site has the
     * summation: D_{ij} x S_{ij}
     *
     * where j is the index of a NN neighbour and D_{ij} = D * r_{ij}, r_{ij}
     * being the position vector connectin the site i with the site j
     *
     * In the continuum, the field can be written as : - nabla X M , which can
     * be discretised with an expression similar to the atomic model one when
     * doing the finite difference approximation of the derivative, thus we
     * only need the DMI vector as:  -D * r_{ij}
     *
     * So, in the loop through neighbours we compute:
     *
     * neighbour          field sum                this gives
     *   -x:      (D / dx) * (+x  X  M)   --> Components in y, z
     *   +x:      (D / dx) * (-x  X  M)   --> %
     *   -y:      (D / dy) * (+y  X  M)   --> Components in x, z
     *   +y:      (D / dy) * (-y  X  M)   --> %
     *   -z:      (D / dz) * (+z  X  M)   --> Components in x, y
     *   +z:      (D / dz) * (-z  X  M)   --> %
     *
     * which gives:
     *      field_x = (D / dx) * (mz[-x] - mz[+x])
     *                  + (D / dz) * (mx[-z] - mx[+x])
     *      ...
     *
     *  which are the first order derivatives.
     *
     *  To compute the DMI for other point groups we will have to find
     *  the vectors and how to discretise the derivative
     *
     */

    double sign;
    double * dmivector = malloc(3 * sizeof(double));
    double DMIc;

    /* Here we iterate through every mesh node */
	#pragma omp parallel for
	for (int i = 0; i < nxyz; i++) {
	    double fx = 0, fy = 0, fz = 0;
	    int idnm = 0;     // Index for the magnetisation matrix
	    int idn = 6 * i; // index for the neighbours

        /* Set a zero field for sites without magnetic material */
	    if (Ms_inv[i] == 0.0){
	        field[3 * i] = 0;
	        field[3 * i + 1] = 0;
	        field[3 * i + 2] = 0;
            energy[i] = 0;
	        continue;
	    }

        /* Here we iterate through the neighbours. Remember:
         * j = 0, 1, 2, 3, 4, 5  --> -x, +x, -y, +y, -z, +z */
        for (int j = 0; j < 6; j++) {

            /* Remember that index=-1 is for sites without material, so
             * we skip those sites */
	        if (ngbs[idn + j] >= 0) {

                /* Since we are iterating through the 6 NN neighbours, we
                 * set the sign for the DMI vector since, for example,
                 * in the x directions, r_{ij} = (+-1, 0, 0)
                 * So, for j=0, 2, 4 (i.e. -x, -y, -z) we use a negatiev sign
                 */
                if (j == 0 || j == 2 || j == 4){ sign = -1; }
                else{ sign = 1; }

                /* Now we set the DMI vectors according to the neighbour position */
                if (j == 0 || j == 1) {
                    DMIc = -D[i] / dx;
                    dmivector[0] = sign, dmivector[1] = 0, dmivector[2] = 0;
                }
                else if (j == 2 || j == 3) {
                    DMIc = -D[i] / dy;
                    dmivector[0] = 0, dmivector[1] = sign, dmivector[2] = 0;
                }
                else if (j == 4 || j == 5) {
                    DMIc = -D[i] / dz;
                    dmivector[0] = 0, dmivector[1] = 0, dmivector[2] = sign;
                }

                /* Magnetisation array index of the neighbouring spin
                 * since ngbs gives the neighbour's index */
	            idnm = 3 * ngbs[idn + j];

                /* Check that the magnetisation of the neighbouring spin
                 * is larger than zero */
                if (Ms_inv[ngbs[idn + j]] > 0){

                    /* We do here:  (D / dx_i) * ( r_{ij} X M_{j} )
                     * The cross_i function gives the i component of the
                     * cross product
                     */

                    /* The x component of the cross product of +-x
                     * times anything is zero (similar for the other comps */
                    if (j != 0 && j != 1) {
                        fx += DMIc * cross_x(dmivector[0], dmivector[1], dmivector[2],
                                             m[idnm], m[idnm + 1], m[idnm + 2]);
                    }
                    if (j != 2 && j != 3) {
                        fy += DMIc * cross_y(dmivector[0], dmivector[1], dmivector[2],
                                             m[idnm], m[idnm + 1], m[idnm + 2]);
                    }
                    if (j != 4 && j != 5) {
                        fz += DMIc * cross_z(dmivector[0], dmivector[1], dmivector[2],
                                             m[idnm], m[idnm + 1], m[idnm + 2]);
                    }
                }
            }
        }

        /* Energy as: (-mu0 * Ms / 2) * [ H_dmi * m ]   */
        energy[i] = -0.5 * (fx * m[3 * i] + fy * m[3 * i + 1]
                            + fz * m[3 * i + 2]);

        /* Update the field H_dmi which has the same structure than *m */
        field[3 * i]     = fx * Ms_inv[i] * MU0_INV;
        field[3 * i + 1] = fy * Ms_inv[i] * MU0_INV;
        field[3 * i + 2] = fz * Ms_inv[i] * MU0_INV;
    }
    free(dmivector);
}


void dmi_field_interfacial(double *m, double *field, double *energy, double *Ms_inv,
                    double *D, double dx, double dy, double dz,
                    int nxyz, int *ngbs) {

    /* In the atomic model, the effective field in the i-th spin site has the
     * summation: D_{ij} x S_{ij}
     *
     * For the Interfacial DMI, D_{ij} has the structure: D_{ij} = r_{ij} X z
     *
     * so the Dzyaloshinskii vectors are IN plane. This function only works
     * along the XY plane since we assume there is a non magnetic material below
     * with a different SOC which gives the DMI for neighbouring in plane spins
     *
     * The computation is similar than before, only that we no longer have
     * a z component. In the continuum the H_dmi field is:
     *          H_dmi = (2 * D / (Ms a a_z)) * ( x X dM/dy  - y X dM/dx )
     * where "a" is the lattice spacing in the plane and "a_z" the system
     * thickness along z
     * This derivative can be obtained doing the cross product
     * as in the atomic model, hence when going through the neighbours
     * loop, we can obtain the field doing the following cross products:
     *
     *   neighbour          field sum
     *     -x:      (D / dx) * (+y  X  M)
     *     +x:      (D / dx) * (-y  X  M)
     *     -y:      (D / dy) * (-x  X  M)
     *     +y       (D / dy) * (+x  X  M)
     * 
     * So, our Dzyaloshinskii vectors can be depicted in a square lattice as
     *
     *                     o  +y
     *                     |
     *                     --> D
     *                ^    | 
     *       -x  o __ | __ o __  | _ o  +x
     *                    |      v
     *                    <--
     *                    |
     *                    o  -y
     *
     * If we start with this picture in the atomic model, we can get the
     * continuum expression when doing the limit  a_x a_y a_z  --> 0
     *
     */

    double sign;
    double * dmivector = malloc(3 * sizeof(double));
    double DMIc;

    /* Here we iterate through every mesh node */
	#pragma omp parallel for
	for (int i = 0; i < nxyz; i++) {
	    double fx = 0, fy = 0, fz = 0;
	    int idnm = 0;     // Index for the magnetisation matrix
	    int idn = 6 * i; // index for the neighbours

        /* Set a zero field for sites without magnetic material */
	    if (Ms_inv[i] == 0.0){
	        field[3 * i] = 0;
	        field[3 * i + 1] = 0;
	        field[3 * i + 2] = 0;
            energy[i] = 0;
	        continue;
	    }

        /* Here we iterate through the neighbours in the XY plane */
        for (int j = 0; j < 4; j++) {

            /* Remember that index=-1 is for sites without material */
	        if (ngbs[idn + j] >= 0) {
                
                /* Since we are iterating through 4 NN neighbours, we
                 * set the sign for the DMI vector since, for example,
                 * in the x directions, (r_{ij} X z) = (0, +-1, 0)
                 * So, for j=1, 2 (i.e. x, -y) we use a negatiev sign
                 * (see the DMI vector details in the above documentation)
                 */
                if (j == 1 || j == 2){ sign = -1; }
                else{ sign = 1; }

                /* Now we set the vectors with the mesh spacing below only
                 * for in plane neighbours */
                if (j == 0 || j == 1) {
                    DMIc = D[i] / dx;
                    dmivector[0] = 0, dmivector[1] = sign, dmivector[2] = 0;
                }
                else if (j == 2 || j == 3) {
                    DMIc = D[i] / dy;
                    dmivector[0] = sign, dmivector[1] = 0, dmivector[2] = 0;
                }

                /* Magnetisation array index of the neighbouring spin
                 * since ngbs gives the neighbour's index */
	            idnm = 3 * ngbs[idn + j];

                /* Check that the magnetisation of the neighbouring spin
                 * is larger than zero */
                if (Ms_inv[ngbs[idn + j]] > 0){

                    fx += DMIc * cross_x(dmivector[0], dmivector[1], dmivector[2],
                                         m[idnm], m[idnm + 1], m[idnm + 2]);
                    fy += DMIc * cross_y(dmivector[0], dmivector[1], dmivector[2],
                                         m[idnm], m[idnm + 1], m[idnm + 2]);
                    fz += DMIc * cross_z(dmivector[0], dmivector[1], dmivector[2],
                                         m[idnm], m[idnm + 1], m[idnm + 2]);
                }
            }
        }

        /* Energy as: (-mu0 * Ms / 2) * [ H_dmi * m ]   */
        energy[i] = -0.5 * (fx * m[3 * i] + fy * m[3 * i + 1]
                            + fz * m[3 * i + 2]);

        /* Update the field H_dmi which has the same structure than *m */
        field[3 * i]     = fx * Ms_inv[i] * MU0_INV;
        field[3 * i + 1] = fy * Ms_inv[i] * MU0_INV;
        field[3 * i + 2] = fz * Ms_inv[i] * MU0_INV;
    }
    free(dmivector);
}
