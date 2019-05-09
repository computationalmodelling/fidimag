#include "clib.h"

/*
 compute the effective exchange field at site i

 H_i = J \sum_<i,j> S_j

 with Hamiltonian

 Hamiltonian = - J \sum_<i,j> S_i \cdot S_j

 Note that the pair <i,j> only run once for each pair.

- NEIGHBOURS in the arguments:

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
field contribution: Sum_j J_j S_j of the neighbours to the exchange interaction

The ngbs array also gives the correct indexes for the spins at periodic
boundaries

 */
void compute_exch_field(double * spin, double * field,
                        double * mu_s_inv,
                        double * energy,
						double Jx, double Jy, double Jz,
                        int * ngbs, int n, int n_ngbs) {

    #pragma omp parallel for
	for (int i = 0; i < n; i++) {

		int id = 0;
		int id_nn = n_ngbs * i; // index for the neighbours

		double fx = 0, fy = 0, fz = 0;

        for (int j = 0; j < 6; j++) {

            if (ngbs[id_nn + j] >= 0) {

                id = 3 * ngbs[id_nn + j];
                fx += Jx * spin[id];
                fy += Jy * spin[id + 1];
                fz += Jz * spin[id + 2];
            }

        }

        field[3 * i] = fx;
        field[3 * i + 1] = fy;
        field[3 * i + 2] = fz;
        energy[i] = -0.5 * (fx * spin[3 * i] + fy * spin[3 * i + 1] +
                            fz * spin[3 * i + 2]);

        // Scale the field to 1/mu_s
        field[3 * i]     *= mu_s_inv[i];
        field[3 * i + 1] *= mu_s_inv[i];
        field[3 * i + 2] *= mu_s_inv[i];
    }
}



double compute_exch_energy(double * spin,
                           double Jx,  double Jy, double Jz,
                           int nx, int ny, int nz, int xperiodic, int yperiodic) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j, k;
	int index, id;
	double Sx, Sy, Sz;
	double energy = 0;

	for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                index = nyz * i + nz * j + k;
                Sx = spin[index];
                Sy = spin[index + n1];
                Sz = spin[index + n2];

                if (i < nx - 1 || xperiodic) {
                    id = index + nyz;
                    if (i == nx-1){
                        id -= n1;
                    }
                    energy += Jx * Sx * spin[id];
                    energy += Jy * Sy * spin[id + n1];
                    energy += Jz * Sz * spin[id + n2];
                }

                if (j < ny - 1 || yperiodic) {
                    id = index + nz;
                    if (j == ny-1){
                        id -= nyz;
                    }
                    energy += Jx * Sx * spin[id];
                    energy += Jy * Sy * spin[id + n1];
                    energy += Jz * Sz * spin[id + n2];
                }

                if (k < nz - 1) {
                    id = index + 1;
                    energy += Jx * Sx * spin[id];
                    energy += Jy * Sy * spin[id + n1];
                    energy += Jz * Sz * spin[id + n2];
                }

            }
        }
	}

	energy = -energy;

	return energy;

}


/*
 compute the effective exchange field at site i

 H_i =  \sum_<i,j> J_{ij} S_j

 with Hamiltonian

 Hamiltonian = - \sum_<i,j> J_{ij} S_i \cdot S_j

 Note that the pair <i,j> only run once for each pair.

 */
void compute_exch_field_spatial(double * spin, double * field,
                                double * mu_s_inv,
                                double * energy,
                                double * J, int * ngbs,
                                int n, int n_ngbs) {

    #pragma omp parallel for
	for (int i = 0; i < n; i++) {

		int id = 0;
		int id_nn = n_ngbs * i; // index for the neighbours

		double fx = 0, fy = 0, fz = 0;

        for (int j = 0; j < 6; j++) {

            int p = id_nn + j;

            if (ngbs[p] >= 0) {

                id = 3 * ngbs[p];
                fx += J[p] * spin[id];
                fy += J[p] * spin[id + 1];
                fz += J[p] * spin[id + 2];
            }

        }

        field[3 * i] = fx;
        field[3 * i + 1] = fy;
        field[3 * i + 2] = fz;
        energy[i] = -0.5 * (fx * spin[3 * i] + fy * spin[3 * i + 1] +
                            fz * spin[3 * i + 2]);

        // Scale the field to 1/mu_s
        field[3 * i]     *= mu_s_inv[i];
        field[3 * i + 1] *= mu_s_inv[i];
        field[3 * i + 2] *= mu_s_inv[i];
    }
}

/* Calculation of Exchange field for up to 8 shells of neighbours
 *
 * J                :: Array with 9 elements: an exchange constant per shell
 *                     Calculation is only up to n_shells (the rest of the elements
 *                     are not used.
 * ngbs             :: array with the neighbours
 * n_ngbs           :: number of neighbours per lattice site
 * n_shells         :: number of specified shells of neighbours
 * n_ngbs_shell     :: number of neighbours per shell. First element is zero.
 *                     For a hex lattice, this array is: [0, 6, 6, 12, ... ]
 * sum_ngbs_shell   :: sum of number of neighbours up to the i-th shell to
 *                     locate the column position for the ngbs of a specific shell.
 *                     First element is zero.
 *                     For a hex lattice, this array is: [0, 6, 12, 24, ...]
 *                     Thus, we can locate ngbs from cols 0-5, 6-11, 12, 23, ... etc
 *
 */
void compute_full_exch_field(double * spin, double * field,
                             double * mu_s_inv,
                             double * energy,
					      	 double J[9], int *ngbs, int n, int n_ngbs,
                             int n_shells, int * n_ngbs_shell,
                             int * sum_ngbs_shell
                             ) {

    #pragma omp parallel for
	for (int i = 0; i < n; i++) {

		int id = 0;
		int id_ngbs = n_ngbs * i; // index for the starting point of neighbours

		double fx = 0, fy = 0, fz = 0;

        for (int sh = 1; sh < n_shells + 1; sh++) {
            for (int j = sum_ngbs_shell[sh - 1]; j < sum_ngbs_shell[sh]; j++) {

                if (ngbs[id_ngbs + j] >= 0) {

                    id = 3 * ngbs[id_ngbs + j];
                    fx += J[sh - 1] * spin[id];
                    fy += J[sh - 1] * spin[id + 1];
                    fz += J[sh - 1] * spin[id + 2];
                }

            }
        }

        field[3 * i] = fx;
        field[3 * i + 1] = fy;
        field[3 * i + 2] = fz;
        energy[i] = -0.5 * (fx * spin[3 * i] + fy * spin[3 * i + 1] +
                            fz * spin[3 * i + 2]);

        // Scale the field to 1/mu_s
        field[3 * i]     *= mu_s_inv[i];
        field[3 * i + 1] *= mu_s_inv[i];
        field[3 * i + 2] *= mu_s_inv[i];
    }
}
