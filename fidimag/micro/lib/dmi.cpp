#include "micro_clib.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
 * D          :: Array with the DMI constant values (see the description on each
 *               function for a detailed explanation)
 *
 * dx, dy, dz :: Mesh spacings in the corresponding directions
 *
 * n          :: Number of mesh nodes
 *
 * ngbs       :: The array of neighbouring spins, which has (6 * n)
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
 *  interfacial DMI we use the interaction found in materials with
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
 * -----------------------------------------------------------------------------
 * -----------------------------------------------------------------------------
 * DMI
 * -----------------------------------------------------------------------------
 * -----------------------------------------------------------------------------
 *
 * In the atomic model, the effective field in the i-th spin site has the
 * summation: D_{ij} x S_{ij}
 *
 * where j is the index of a NN neighbour and D_{ij} = D * r_{ij}, r_{ij}
 * being the position vector connectin the site i with the site j
 *
 * ------------------------------------------------------------------------
 * BULK DMI
 * ------------------------------------------------------------------------
 *
 * In the continuum, the field can be written as (M = Ms m):
 *                                           ->
 *                  - (2 D / mu0 Ms) nabla X m ,
 *
 * which can be discretised with an expression similar to the atomic model
 * one when doing the finite difference approximation of the derivative,
 * thus we only need the DMI vector as:  -D * r_{ij}. Specifically, we have
 * that:
 *
 *      ---->                  ->                       ^         ^
 *      field = - 2 D  nabla X m = -D [ (d m_z - d m_y) x + (...) y +  ]
 *               ----                    -----   -----
 *               mu0 Ms                   d y     d z
 *
 * Discretising the derivatives, the field components are:
 *
 *      field_x = - (2 D / mu0 Ms) [ (1 / 2 dy) * (mz[+y] - mz[-y]) -
 *                                   (1 / 2 dz) * (my[+z] - my[-z])
 *                                 ]
 *
 *      field_y = - (2 D / mu0 Ms) [ (1 / 2 dz) * (mx[+z] - mx[-z]) -
 *                                   (1 / 2 dx) * (mz[+x] - mz[-x])
 *                                 ]
 *      ...
 *
 * where +x, -x, etc are the neighbours at position +x, -x, etc
 * We now can collect the terms involving the neighbours at
 * +x, -x, +y, -y, +z and -z.
 * For example the field H where the +x neighbour is involved is:
 *
 *  ->                                 ^          ^
 *  H(+x) = (-2 D / mu0 Ms) * ( my(+x) z - mz(+x) y ) * (1 / 2 dx)
 *
 *                                 ^     ->
 *        = (- D / mu0 Ms dx) * (  x  X  m  )
 *
 *  In general, the field for every component will be the cross product of
 *  r_ij with m, where r_ij is the vector towards the j-th neighbour, as in
 *  the discrete case. This term is divided by the mesh discretisation
 *  (dx, dy or dz) in the corresponding direction
 *
 *  The DMI vector norms are given by the *D array. If our simulation has
 *  n mesh nodes, then the D array len is n. The order is the same than the NNs
 *  array, i.e.
 *
 *      D = [D_0, D_1, ...]
 *
 *  where D_i means the DMI vector norm of the i-th mesh node.
 *
 *  NOTE:
 *  To compute the DMI for other point groups we will have to find
 *  the vectors and how to discretise the derivative
 *
 *  MULTIPLE DMI --------------------------------------------------------------
 *
 *  Multiple DMI vectors can be passed in the dmi_vector array. For ex, with
 *  two DMIs:
 *
 *         dmi_vector = [ Dx(-x) Dy(-x) Dz(-x) -- DMI vector 1 per every NN
 *                        Dx(+x) Dy(+x) Dz(+x) 
 *                        Dx(-y) Dy(-y) Dz(-y) 
 *                        ... 
 *                        Dx(+z) Dy(+z) Dz(+z) --
 *                        Dx(-x) Dy(-x) Dz(-x) -- DMI vector 2 per every NN
 *                        Dx(+x) Dy(+x) Dz(+x) 
 *                        Dx(-y) Dy(-y) Dz(-y) 
 *                        ... 
 *                        Dx(+z) Dy(+z) Dz(+z) 
 *                      ]
 *
 * and DMI constants are passed inerspersed in the D array per every node
 *   
 *    D = [ D1_0 D2_0 D1_1 D2_1 ... D1_N D2_N]
 *                      
 *  where D1_i is the 1st DMI constant of the i-mesh node
 * 
 * ------------------------------------------------------------------------
 * INTERFACIAL DMI
 * ------------------------------------------------------------------------
 *
 * For the Interfacial DMI, we use the following energy density structure
 * with Lifshitz invariants (see Rohart et al. Phys. Rev. B 88, 184422))
 *
 *      w_DM = D ( L_{xz}^{(x)} + L_{yz}^{(y)} )
 *
 * Hence, using a variational derivative, the field can be computed as
 *
 *                       /        ->           ->	\
 *      field = - 2  D  |  ^     dm     ^     dm    |
 *                ----  |  y  X  --  -  x  X  --    |
 *               mu0 Ms  \       dx           dy   /
 *
 * Using finite differences for the derivatives, we can get the components
 * contributed to the field by the neighbours in the +x, -x, +y and -y
 * directions. For instance,
 *
 *                        /				\
 *  field(+x)  = -2  D   |  mz(+x) ^      mx(+x) ^  |
 *                ----   |  -----  x  -   -----  z  |
 *               mu0 Ms   \  2 dx          2 dx    /
 *
 *               - D     1     ^   ->
 *             =  ----   --  ( y X M(+x) )
 *               mu0 Ms  dx
 *
 *
 * In general, the derivative will be given as in the discrete spin model
 * by using the Dzyaloshinskii vector as
 *
 *                 ^   ->
 *          D_ij = z X r_ij
 *
 * and dividing by the mesh discretisation in the direction of the
 * neighbour.
 * This DMI vector convention is given by [Yang et al. PRL 115, 267210]
 *
 * The Dzyaloshinskii vectors are IN plane. This function only works along
 * the XY plane since we assume there is a non magnetic material below with
 * a different SOC which gives the DMI for neighbouring in plane spins
 *
 * As in the atomic model, when going through the neighbours
 * loop, we can obtain the field doing the following cross products:
 *
 *   neighbour          field sum
 *     -x:      (D / dx) * (-y  X  M)
 *     +x:      (D / dx) * (+y  X  M)
 *     -y:      (D / dy) * (+x  X  M)
 *     +y       (D / dy) * (-x  X  M)
 *     +z       0
 *     -z       0
 *
 * So, our Dzyaloshinskii vectors can be depicted in the cuboid mesh as
 *
 *                     o  +y
 *                     |
 *                   <--  D
 *                     |     ^
 *       -x  o __ | __ o __  | _ o  +x
 *                v    |
 *                    -->
 *                     |
 *                     o  -y
 *
 * If we start with this picture in the atomic model, we can get the
 * continuum expression when doing the limit  a_x a_y a_z  --> 0
 */


// ----------------------------------------------------------------------------


/* A function to compute any DMI given the corresponding DM vectors
 * obtained from discretising the Lifshitz invariants using the
 * finite differences:
 *
 * The calculation of the field for the neighbour in the i-direction is done as:
 *
 *      H_DMI = -2 D     1       ->        ->
 *              ----    ----     x_{i}  X  m
 *              mu0 Ms  2 dxi
 *
 * where x_{i} is the DM vector orientation in the i-direction and dxi
 * is the mesh spacing in that direction, i.e. for a neighbour in the
 * +y direction it is dy
 * These orientations are specified in the dmi_vector array
 *
 * DM vectors for bulk and interfacial DMI are explained at the beginning
 * of this library
 *
 * dmi_vector   :: an array with a vector for every nearest neighbour from
 *                 the finite differences model. For DMIs defined in 2D
 *                 the last 6 terms are set to zero
 *
*/
void dmi_field(double * m, double * field,
               double * energy, double * Ms_inv,
               double * D, int n_dmis, 
               double *dmi_vector,
               double dx, double dy, double dz,
               int n, int * ngbs) {

    /* These are for the DMI prefactor or coefficient */
    double dxs[6] = {dx, dx, dy, dy, dz, dz};

    /* Here we iterate through every mesh node */
    for (int i = 0; i < n; i++) {
        double DMIc;
        double fx = 0, fy = 0, fz = 0;
        int idnm;     // Index for the magnetisation matrix
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
         * j = 0, 1, 2, 3, 4, 5  --> -x, +x, -y, +y, -z, +z
         Previously we had n_dmi_neighbours as a parameter,
         but now we just skip if the vector entries are zero
         instead, because we can have cases such as in D_n
         where the D1 component has -x,+x,-y,+y,0,0
         and the D2 component has 0,0,0,0,-z,+z.
        */
        for (int j = 0; j < 6; j++) {
            /* Add the DMI field x times for every DMI constant */
            for (int k = 0; k < n_dmis; k++) {
                
                // starting index of the DMI vector for this neighbour (j)
                // (remember we have 18 comps of dmi_vector per DMI constant)
                int ngbr_idx_D = k * 18 + 3 * j;

                /* We skip the neighbour if
                   (a) it doesn't exist (ngbs[idn + j] = -1)
                   (b) there is no material there
                   (c) DMI value is zero there
                */
                if ((ngbs[idn + j] != -1) && (Ms_inv[ngbs[idn + j]] !=  0)) {
                /* We do here:  -(D / dx_i) * ( r_{ij} X M_{j} )
                 * The cross_i function gives the i component of the
                 * cross product. The coefficient is computed according
                 * to the DMI strength of the current lattice site.
                 * For the denominator, for example, if j=2 or 3, then
                 * dxs[j] = dy
                 */

                    /* check the DMI vector exists */
                    if (dmi_vector[ngbr_idx_D    ] != 0 ||
                        dmi_vector[ngbr_idx_D + 1] != 0 ||
                        dmi_vector[ngbr_idx_D + 2] != 0   ) {
                        /* Notice the x component of the cross product of +-x
                         * times anything is zero (similar for the other comps) */

                        // Get DMI constant from neighbour and scale
                        // DMI components are in consecutive order per neighbour,
                        // i.e. if 2 DMI consts: [D0 D1 D0 D1 .... ]
                        DMIc = -D[n_dmis * ngbs[idn + j] + k] / dxs[j];

                        idnm = 3 * ngbs[idn + j]; // index for magnetisation

                        fx += DMIc * cross_x(dmi_vector[ngbr_idx_D],
                                             dmi_vector[ngbr_idx_D + 1],
                                             dmi_vector[ngbr_idx_D + 2],
                                             m[idnm], m[idnm + 1], m[idnm + 2]);
                        fy += DMIc * cross_y(dmi_vector[ngbr_idx_D],
                                             dmi_vector[ngbr_idx_D + 1],
                                             dmi_vector[ngbr_idx_D + 2],
                                             m[idnm], m[idnm + 1], m[idnm + 2]);
                        fz += DMIc * cross_z(dmi_vector[ngbr_idx_D],
                                             dmi_vector[ngbr_idx_D + 1],
                                             dmi_vector[ngbr_idx_D + 2],
                                             m[idnm], m[idnm + 1], m[idnm + 2]);

                    }
                }
            }  // Close for loop through n of DMI constants
        }  // Close for loop through neighbours per mesh site

        /* Energy as: (-mu0 * Ms / 2) * [ H_dmi * m ]   */
        energy[i] = -0.5 * (fx * m[3 * i] + fy * m[3 * i + 1]
        		+ fz * m[3 * i + 2]);

        /* Update the field H_dmi which has the same structure than *m */
        field[3 * i]     = fx * Ms_inv[i] * MU0_INV;
        field[3 * i + 1] = fy * Ms_inv[i] * MU0_INV;
        field[3 * i + 2] = fz * Ms_inv[i] * MU0_INV;
    }
}
