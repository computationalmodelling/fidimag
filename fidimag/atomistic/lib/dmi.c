#include "clib.h"
#include "math.h"
#include "stdlib.h"

inline void normalised_cross_product(double *a, int ai, double *b, int bi,
                                     double *c, int ci) {
    /* Compute the cross product: a X b
     * and generate the result in the array c
     *
     * using three array components starting at
     * the ai, bi or ci-th position for a, b and c respectively,
     * e.g. a[ai], a[ai + 1], a[ai + 2],
     * for the x, y and z components
     */
    double norm;

    c[ci]     = a[ai + 1] * b[bi + 2] - a[ai + 2] * b[bi + 1];
    c[ci + 1] = a[ai + 2] * b[bi]     - a[ai]     * b[bi + 2];
    c[ci + 2] = a[ai]     * b[bi + 1] - a[ai + 1] * b[bi];

    norm = sqrt(c[ci] * c[ci] + c[ci + 1] * c[ci + 1] + c[ci + 2] * c[ci + 2]);

    if (abs(norm) < 1 - 1e-6 || abs(norm) > 1 + 1e-6){
        /* We will assume the norm is not zero */
        c[ci] = c[ci] / norm;
        c[ci + 1] = c[ci + 1] / norm;
        c[ci + 2] = c[ci + 2] / norm;
    }
}

void dmi_field_bulk(double *spin, double *field,
                    double *energy, double D, int *ngbs, int nxyz) {

    /* Bulk DMI field and energy computation
     *
     * ngbs[] contains the *indexes* of the neighbours
     * in the following order:
     *      -x, +x, -y, +y, -z, +z
     *
     * for every spin. It is -1 for boundaries.
     * The array is like:
     *      | 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, 1-x, 1+x, 1-y, ...  |
     *        i=0                           i=1                ...
     *
     * where  0-y  is the index of the neighbour of the 0th spin,
     * in the -y direction, for example
     *
     * Thus, for every neighbour ( ngbs[i + j], j=0,1,...5 )
     * we compute the field contribution.
     *
     * The value of ngbs[] at sites with Ms = 0, is negative
     *
     * The ngbs array also gives the correct indexes for the spins
     * at periodic boundaries
     *
     * The Bulk DMI Dzyaloshinskii vectors are in the same direction than the
     * r_ij vectors which connect two neighbouring sites, pointing from the
     * lattice site *i* towards *j*
     *
     * The bulk DMI is defined in 3 dimensions
     *
     * Then the DMI field is computed as :
     *
     *      Sum_j ( D_ij X S_j ) = Sum_j ( r_ij X S_j )
     *
     * for every spin *i*
     *
     * --- Check this: can we use a SC lattice for this bulk DMI expression?
     *     Since, for example, MnSi crystal structure is more complex 
     */

	#pragma omp parallel for
	for (int i = 0; i < nxyz; i++) {

		int id = 0;
		int idv = 6 * i; // index for the neighbours

		double fx = 0, fy = 0, fz = 0;

		if (ngbs[idv]>=0) { // neighbour at x-1
			id = 3*ngbs[idv];
			fx += D*cross_x(-1,0,0,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(-1,0,0,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(-1,0,0,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+1]>=0) { // neighbour x+1
			id = 3*ngbs[idv+1];
			fx += D*cross_x(1,0,0,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(1,0,0,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(1,0,0,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+2]>=0) { // neighbour at y-1
			id = 3*ngbs[idv+2];
			fx += D*cross_x(0,-1,0,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(0,-1,0,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(0,-1,0,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+3]>=0) { // neighbour at y+1
			id = 3*ngbs[idv+3];
			fx += D*cross_x(0,1,0,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(0,1,0,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(0,1,0,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+4]>=0) { // neighbour at z-1
			id = 3*ngbs[idv+4];
			fx += D*cross_x(0,0,-1,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(0,0,-1,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(0,0,-1,spin[id],spin[id+1],spin[id+2]);
		}

		if (ngbs[idv+5]>=0) { // neighbour at z+1
			id = 3*ngbs[idv+5];
			fx += D*cross_x(0,0,1,spin[id],spin[id+1],spin[id+2]);
			fy += D*cross_y(0,0,1,spin[id],spin[id+1],spin[id+2]);
			fz += D*cross_z(0,0,1,spin[id],spin[id+1],spin[id+2]);
		}

		field[3 * i] = fx;
        field[3 * i + 1] = fy;
        field[3 * i + 2] = fz;

        energy[i] = -0.5 * (fx * spin[3 * i] +
                            fy * spin[3 * i + 1] +
                            fz * spin[3 * i + 2]
                            );
    }
}

void dmi_field_interfacial_atomistic(double *spin,
    double *field, double *energy,
    double D, int *ngbs, int nxyz, int nneighbours, double *r, int rdim) {

    /* Interfacial DMI field and energy computation
     *
     * ngbs[] contains the *indexes* of the neighbours
     * for different lattices. The number of neighbours is specified with the
     * nneighbours integer
     *
     * For example, in a simple cubic (SC) lattice,
     * the array is like:
     *      | 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, 1-x, 1+x, 1-y, ...  |
     *        i=0                           i=1                ...
     *
     * where  0-y  is the index of the neighbour of the 0th spin,
     * in the -y direction, for example
     *
     * Thus, for every neighbour of the i-th spin ( ngbs[i + j], j=0,1,...5 )
     * we compute the field contribution
     * The value of ngbs[] at sites with Ms = 0, is negative
     *
     * The ngbs array also gives the correct indexes for the spins
     * at periodic boundaries
     *
     * The Interfacial DMI needs the direction of the Dzyaloshinskii vectors,
     * which are obtained as
     *
     *      D_ij = r_ij  X  +z
     *
     * where r_ij is the vector connecting two neighbouring sites.
     * The r_ij vectors will be computed using the position vectors stored
     * in the *r array, which is rdim-dimensional (hexagonal lattice is
     * 2D and SC lattice is 3D)
     * 
     * For instance, in a SC 2D lattice:
     * [ lattice site X,  neighbours O ]
     *
     *                    O
     *                    .
     *                   --->  D_ij
     *                    .
     *              ^     .                  ^ y
     *          O . | . . X . . | . O        |
     *                    .     v            |---> x
     *                    .
     *                   <---
     *                    .
     *                    O
     *
     * (this DMI is only defined in two dimensions--> interfaces, thin films)
     * Then the DMI field is computed as :  Sum_j ( D_ij X S_j )
     * for every spin *i*
     *
     * In a hexagonal lattice we have 6 neighbours
     *
     * We assume a constant DMI vector magnitude
     */


  	#pragma omp parallel for
	for (int i = 0; i < nxyz; i++) {

        /* Generate the direction of the Dzyaloshinskii vectors for the
         * 6 nearest neighbours
         * Since we have 3 directions for every vector, we have in total
         * 18 entries
         */
        double * DMI_vec = malloc(6 * 3 * sizeof(double));
        for(int j; j < 6; i++) DMI_vec[j] = 0;
        
        /* Normalised z vector */
        double * z_vec = malloc(3 * sizeof(double));
        z_vec[0] = 0, z_vec[1] = 0, z_vec[2] = 1;
        
        /* Array to compute the position vector fromt the i-th site 
         * to the j-th site */
        double * rij = malloc(3 * sizeof(double));

		int idn = 6 * i; // index for the NNs
		int idnm = 0; // index for the magnetisation components of the NNs

		double fx = 0, fy = 0, fz = 0;

        /* Now we compute for every neighbour */
        for(int j = 0; j < nneighbours; j++){
            /* Check that Ms != 0 */
            if (ngbs[idn + j] >= 0) {

                /* We are assuming the Dzyaloshinskii vectors are the same in
                 * equivalent directions for every lattice site, so we only
                 * update their values when necessary  */
                if (DMI_vec[j] != 0 && DMI_vec[j + 1] != 0 && DMI_vec[j + 2] != 0){
                    /* Position vectors. Since we assume a 2D lattice, we set 
                     * the z component as zero */
                    rij[0] = r[ngbs[idn + j]]     - r[rdim * i];
                    rij[1] = r[ngbs[idn + j] + 1] - r[rdim * i + 1];
                    /* rij[2] = r[ngbs[idn + j] + 2] - r[rdim * i + 2]; */
                    rij[2] = 0;

                    /* Now we compute: r_ij X z We save the results in the j-th
                     * neighbour of the DMI vectors array */
                    normalised_cross_product(rij, 0, z_vec, 0, 
                                             DMI_vec, 3 * j);
                }

                /* Now we add the field contribution of the j-th
                 * neighbour for the i-th spin: D_ij X S_j */
                idnm = 3 * ngbs[idn];
                fx += D * cross_x(DMI_vec[3 * j], DMI_vec[3 * j + 1], DMI_vec[3 * j + 2],
                                  spin[idnm], spin[idnm + 1], spin[idnm + 2]);
                fy += D * cross_y(DMI_vec[3 * j], DMI_vec[3 * j + 1], DMI_vec[3 * j + 2],
                                  spin[idnm], spin[idnm + 1], spin[idnm + 2]);
                fz += D * cross_z(DMI_vec[3 * j], DMI_vec[3 * j + 1], DMI_vec[3 * j + 2],
                                  spin[idnm], spin[idnm + 1], spin[idnm + 2]);
            }
        }

        field[3 * i] = fx;
  	    field[3 * i + 1] = fy;
  	    field[3 * i + 2] = fz;

        // TODO: check whether the energy is correct or not.
        /* Avoid second counting with 1/2 */
  	    energy[i] = -0.5 * (fx * spin[3 * i] +
                            fy * spin[3 * i + 1]+
                            fz * spin[3 * i + 2]
                            );

        free(DMI_vec);
        free(z_vec);
        free(rij);
    }
}


inline double single_energy_x(double D, double Si[3], double Sj[3]){
    double tx = D*(-Si[2]*Sj[1]+Si[1]*Sj[2]);
    return tx;
}

inline double single_energy_y(double D, double Si[3], double Sj[3]){
    double ty = D*(Si[2]*Sj[0]-Si[0]*Sj[2]);
    return ty;
}

inline double single_energy_z(double D, double Si[3], double Sj[3]){
    double tz = D*(-Si[1]*Sj[0]+Si[0]*Sj[1]);
    return tz;
}

double dmi_energy(double *spin, double D, int nx, int ny, int nz, int xperiodic, int yperiodic) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int index, id;

	double energy = 0;

    double S_i[3],S_j[3];

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				index = nyz * i + nz * j + k;
				S_i[0] = spin[index];
				S_i[1] = spin[index + n1];
				S_i[2] = spin[index + n2];

				if (i < nx - 1 || xperiodic) {
					id = index + nyz;
                    if (i == nx -1){
                        id -= n1;
                    }
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_x(D,S_i,S_j);
				}

				if (j < ny - 1 || yperiodic) {
					id = index + nz;
                    if (j == ny-1){
                        id -= nyz;
                    }
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_y(D,S_i,S_j);
				}

				if (k < nz - 1) {
					id = index + 1;
                    S_j[0] = spin[id];
                    S_j[1] = spin[id + n1];
                    S_j[2] = spin[id + n2];
					energy += single_energy_z(D,S_i,S_j);
				}

			}
		}
	}

	return energy;

}
