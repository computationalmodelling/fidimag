#include "clib.h"
#include "math.h"
#include "stdlib.h"

/* NOTES about the neighbours array in the arguments:

The first 6 elements of the ngbs[] are the indexes of the nearest neighbours in
the following order:

     -x, +x, -y, +y, -z, +z

for every spin. It is -1 for boundaries.  The array is like:

                                      __here are next neighbour indexes
                                      | 
     | 0-x, 0+x, 0-y, 0+y, 0-z, 0+z, ... 1-x, 1+x, 1-y, ...  |
       i=0                               i=1                ...

where  0-y  is the index of the neighbour of the 0th spin, in the -y direction,
for example

Thus, for every nearest neighbour ( ngbs[i + j], j=0,1,...5 ) we compute the
field contribution from the corresponding DMI

The ngbs array also gives the correct indexes for the spins at periodic
boundaries

*/

void dmi_field_bulk(double *spin, double *field,
                    double *energy, double *_D, int *ngbs, int nxyz, int n_ngbs) {

    /* Bulk DMI field and energy computation
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
     * NOTES:
     *
     * 1. Check this: can we use a SC lattice for this bulk DMI expression?
     *    Since, for example, MnSi crystal structure is more complex
     *
     * 2. This function is not defined for a hexagonal lattice, we can
     *    generalise the code to solve this
     *      
     */


    /* The DMI vectors are the same according to the neighbours
     * positions. Thus we set them here to avoid compute them
     * every time in the loop . So, if we want the j-th NN,
     * the DMI vector will be 
     * (D_x, D_y, D_z) --> ( dmivector[3 * j] ,
     *                       dmivector[3 * j + 1],
     *                       dmivector[3 * j + 2] )
     */
    double dmivector[18] = {-1,  0,  0,
                             1,  0,  0,
                             0, -1,  0,
                             0,  1,  0,
                             0,  0, -1,
                             0,  0,  1
                             };

	#pragma omp parallel for shared(dmivector)
	for (int i = 0; i < nxyz; i++) {

		int id = 0;
		int id_nn = n_ngbs * i; // index for the neighbours

		double fx = 0, fy = 0, fz = 0;
        double D=0;

        // Compute the DMI contribution for every NEAREST neighbour:
        // -x, +x, -y, +y, -z, +z
        for (int j = 0; j < 6; j++){
            if (ngbs[id_nn + j] >= 0) {
                id = 3 * ngbs[id_nn + j];
                D = _D[id_nn + j];
                fx += D * cross_x(dmivector[3 * j], 
                                  dmivector[3 * j + 1], 
                                  dmivector[3 * j + 2], 
                                  spin[id], spin[id+1], spin[id+2]);
                fy += D * cross_y(dmivector[3 * j], 
                                  dmivector[3 * j + 1], 
                                  dmivector[3 * j + 2], 
                                  spin[id], spin[id+1], spin[id+2]);
                fz += D * cross_z(dmivector[3 * j], 
                                  dmivector[3 * j + 1], 
                                  dmivector[3 * j + 2], 
                                  spin[id], spin[id+1], spin[id+2]);
            }
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

void dmi_field_interfacial_atomistic(double *spin, double *field, double *energy,
    double D, int *ngbs, int n, int n_ngbs, int n_ngbs_dmi, double *DMI_vec) {
    
    /* Interfacial DMI field and energy computation
     *
     * ngbs[] contains the *indexes* of the neighbours
     * for different lattices. The number of neighbours is specified with the
     * n_ngbs integer.
     * The number of neighbours used to compute the 2D interfacial DMI is
     * given by the n_dmi_ngbs (e.g. 4 in a square lattice, 6 in hex lattice)
     *
     * The Interfacial DMI needs the direction of the Dzyaloshinskii vectors,
     * which are obtained as
     *
     *      D_ij = r_ij  X  +z
     *
     * where r_ij is the vector connecting two neighbouring sites.
     * This vectors are computed in Python and passed here as the DMI_vec array,
     * whose entries are according to the neighbours matrix order, i.e.
     * the (DMI_vec[3 * j], DMI_vec[3 * j + 1], DMI_vec[3 * j + 2]) vector
     * are the Dij vector components between the i-th site and the j-th neighbour
     * (we assume the vectors are the same between NNs for every lattice site)
     *
     * 
     * For instance, in a SC 2D lattice, the X -th spin has the DMI vectors:
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
     * so DMI_vec = {0, 1, 0, 0, -1, 0, -1, 0, 0, 1, 0, 0}
     * NN:          j=0      j=1        j=2      j=3
     *
     * This DMI is only defined in two dimensions--> interfaces, thin films
     * Then the DMI field is computed as :  Sum_j ( D_ij X S_j )
     * for every spin *i*
     *
     * In a hexagonal lattice we have 6 neighbours
     *
     * We assume a constant DMI vector magnitude
     */
    
  	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
        
		int idn = n_ngbs * i; // index for the NNs
		int idnm = 0;         // index for the magnetisation components of the NNs

		double fx = 0, fy = 0, fz = 0;

        /* Now we compute the field contribution from every NEAREST neighbour in the plane */
        for(int j = 0; j < n_ngbs_dmi; j++){
            /* Check that Ms != 0 */
            if (ngbs[idn + j] >= 0) {

                /* Now we add the field contribution of the j-th
                 * neighbour for the i-th spin: D_ij X S_j */
                idnm = 3 * ngbs[idn + j];  // Magnetisation component index of the j-th NN
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
