#include "common_clib.h"

/* Calculation of the STT field using the Zhang Li formalism, which
 * is
 *          ( \vec{j}_s \cdot \nabla ) \vec{m}
 *
 * where j_s is the current vector in 3D
 * For this computation we use the neighbours matrix to make easier the
 * derivative calculation at the boundaries.
 * ** For now this only works with the CUBOID mesh
 *
 * We have two derivatives, d/dx and d/dy which we approximate using a
 * second order finite difference. Thus, we must distinguish the boundary
 * cases where there is a single neighbour:
 *
 *  - With two neighbours, i.e. in the inner region of the mesh,
 *  the derivatives are simply:
 *          f(x + 1, y, z) - f(x - 1, y, z) / 2 dx
 *          f(x, y + 1, z) - f(x, y - 1, z) / 2 dy
 *          f(x, y, z + 1) - f(x, y, z - 1) / 2 dz
 *
 *  where x+-1 is the neighbour to the left or right.
 *
 *  - With a single neighbour, it  must be changed to
 *          f(x, y, z) - f(x - 1, y, z) / dx    --> no NN to the right
 *    OR
 *          f(x + 1, y, z) - f(x, y, z) / dx    --> no NN to the left
 *                                                  (e.g. no PBC for the 0th
 *                                                   spin in a mesh)
 *
 * In the code, we denote nn_x1 for the index of the NN at the left of the
 * i-th site and nn_x2 for the NN to the right. These variables are simply
 * the i-th index (in the corresponding cases) when there is a single NN
 *
 * Similar for y and z
 *
 * n is the number of spins or lattice sites in the system
 *
 * Vector Fields f have 3 * n entries (e.g. magnetisation=spin, stt field)
 * in the order: [fx0, fy0, fz0, fx1, fy1, fz1, fx2, ...]
 *
 * Scalar fields have n entries
 *
 * The nearest neighbour matrix contains 6 * n indexes, that are the indexes
 * of the 6 NNs for every lattice site.
 * The order of this array is:
 *
 *      ngbs = [ x0-1 x0+1 y0-1 y0+1 z0-1 z0+1 x1-1 y1-1 ...]
 *
 *  where xi-1 is the index of the NN in the -x direction of the i-th site
 *  Neighbouring sites where there is no material, has index -1
 *
 */
void compute_stt_field_c(double *restrict spin, double *restrict field, double *restrict jx, double *restrict jy, double *restrict jz,
		double dx, double dy, double dz, int *restrict ngbs, int n) {

    //#pragma omp parallel for
	for (int i = 0; i < 3 * n; i++) {
		field[i] = 0;
	}

    #pragma omp parallel for
    /* Iterate through every lattice site */
    for (int i = 0; i < n; i++){

        /* Starting index for the NNs of the i-th site
         * i+0, i+1, i+2, i+3 ...  --> -x, +x, -y, +y ...
         */
        int nn_i = 6 * i;
        double factor_x, factor_y, factor_z;
        int nn_x1, nn_x2, nn_y1, nn_y2, nn_z1, nn_z2;

        /* Here we distinguish if there are 2 NNs, no NN in the
         * -x direction, or no NN in the +x direction,
         *  or no NN  at all in the x dir.
         *  In the latest case, make factor_x equal to zero to avoid
         *  summing field to the for loop
         */
        if(ngbs[nn_i] != -1 && ngbs[nn_i + 1] != -1) {
            factor_x = 2;
            nn_x1 = ngbs[nn_i];
            nn_x2 = ngbs[nn_i + 1];
        // Here there is no NN to the right so we make f(x) - f(x-1)
        } else if(ngbs[nn_i + 1] == -1 && ngbs[nn_i] != -1){
            factor_x = 1;
            nn_x1 = ngbs[nn_i];
            nn_x2 = i;
        // Here there is no NN to the left so we make f(x + 1) - f(x)
        } else if(ngbs[nn_i] == -1 && ngbs[nn_i + 1] != -1){
            factor_x = 1;
            nn_x1 = i;
            nn_x2 = ngbs[nn_i + 1];
        } else {
            factor_x = 0;
        }

        /* Here we loop to sum to the x, y, z (3*i+0, 3*i+1, 3*i+2) component
         * of the field for the i-th spin
         * jx is a scalar field, so it only has n entries
         * This calculation is:  jx[i] * d m[i] / dx
         * */
        if (factor_x){
            for(int j = 0; j < 3; j++){
                field[3 * i + j] += jx[i] * (spin[3 * nn_x2 + j]
                                             - spin[3 * nn_x1 + j]) / (factor_x * dx);
            }
        }

        // We do the same along the y direction
        if(ngbs[nn_i + 2] != -1 && ngbs[nn_i + 3] != -1) {
            factor_y = 2;
            nn_y1 = ngbs[nn_i + 2];
            nn_y2 = ngbs[nn_i + 3];
        } else if(ngbs[nn_i + 3] == -1 && ngbs[nn_i + 2] != -1){
            factor_y = 1;
            nn_y1 = ngbs[nn_i + 2];
            nn_y2 = i;
        } else if(ngbs[nn_i + 2] == -1 && ngbs[nn_i + 3] != -1){
            factor_y = 1;
            nn_y1 = i;
            nn_y2 = ngbs[nn_i + 3];
        } else {
            factor_y = 0;
        }

        if (factor_y){
            for(int j = 0; j < 3; j++){
                field[3 * i + j] += jy[i] * (spin[3 * nn_y2 + j]
                                             - spin[3 * nn_y1 + j]) / (factor_y * dy);
            }
        }


        // We do the same along the z direction
        if(ngbs[nn_i + 4] >= 0 && ngbs[nn_i + 5] >= 0) {
            factor_z = 2;
            nn_z1 = ngbs[nn_i + 4];
            nn_z2 = ngbs[nn_i + 5];
        } else if(ngbs[nn_i + 4] >= 0 && ngbs[nn_i + 5] < 0){
            factor_z = 1;
            nn_z1 = ngbs[nn_i + 4];
            nn_z2 = i;
        } else if(ngbs[nn_i + 4] < 0 && ngbs[nn_i + 5] >= 0 ){
            factor_z = 1;
            nn_z1 = i;
            nn_z2 = ngbs[nn_i + 5];
        } else {
            factor_z = 0;
        }

        if (factor_z){
            for(int j = 0; j < 3; j++){
                field[3 * i + j] += jz[i] * (spin[3 * nn_z2 + j]
                                             - spin[3 * nn_z1 + j]) / (factor_z * dz);
            }
        }

    }
}


void llg_stt_rhs(double *restrict dm_dt, double *restrict m, double *restrict h, double *restrict h_stt,
		double *restrict alpha, double beta, double u0, double gamma, int n) {

	#pragma omp parallel for
	for (int index = 0; index < n; index++) {
	    int i = 3 * index;
	    int j = 3 * index + 1;
	    int k = 3 * index + 2;

	    double coeff = -gamma / (1 + alpha[index] * alpha[index]);

	    double mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
	    double mh = m[i] * h[i] + m[j] * h[j] + m[k] * h[k];

            //hp=mm.h-mh.m=-mx(mxh)
            double hpi = mm*h[i] - mh*m[i];
            double hpj = mm*h[j] - mh*m[j];
            double hpk = mm*h[k] - mh*m[k];

	    double mth0 = cross_x(m[i], m[j], m[k], hpi, hpj, hpk);
	    double mth1 = cross_y(m[i], m[j], m[k], hpi, hpj, hpk);
	    double mth2 = cross_z(m[i], m[j], m[k], hpi, hpj, hpk);

	    dm_dt[i] = coeff * (mth0 - hpi * alpha[index]);
	    dm_dt[j] = coeff * (mth1 - hpj * alpha[index]);
	    dm_dt[k] = coeff * (mth2 - hpk * alpha[index]);

	    //the above part is standard LLG equation.

	    double coeff_stt = u0 / (1 + alpha[index] * alpha[index]);

	    double mht = m[i] * h_stt[i] + m[j] * h_stt[j] + m[k] * h_stt[k];

	    hpi = mm*h_stt[i] - mht * m[i];
	    hpj = mm*h_stt[j] - mht * m[j];
	    hpk = mm*h_stt[k] - mht * m[k];

	    mth0 = cross_x(m[i], m[j], m[k], hpi, hpj, hpk);
	    mth1 = cross_y(m[i], m[j], m[k], hpi, hpj, hpk);
	    mth2 = cross_z(m[i], m[j], m[k], hpi, hpj, hpk);

	    dm_dt[i] += coeff_stt * ((1 + alpha[index] * beta) * hpi
				     - (beta - alpha[index]) * mth0);
	    dm_dt[j] += coeff_stt * ((1 + alpha[index] * beta) * hpj
				     - (beta - alpha[index]) * mth1);
	    dm_dt[k] += coeff_stt * ((1 + alpha[index] * beta) * hpk
				     - (beta - alpha[index]) * mth2);

	    double c = 6 * sqrt(dm_dt[i] * dm_dt[i] +
				dm_dt[j] * dm_dt[j] +
				dm_dt[k]* dm_dt[k]);

	    dm_dt[i] += c * (1 - mm) * m[i];
	    dm_dt[j] += c * (1 - mm) * m[j];
	    dm_dt[k] += c * (1 - mm) * m[k];

	}

}


void llg_stt_cpp(double *restrict dm_dt, double *restrict m, double *restrict h, double *restrict p,
		double *restrict alpha, int *restrict pins, double *restrict a_J, double beta, double gamma, int n) {

	#pragma omp parallel for
	for (int index = 0; index < n; index++) {
	    int i = 3 * index;
	    int j = 3 * index + 1;
	    int k = 3 * index + 2;

	    if (pins[index]>0){
		dm_dt[i] = 0;
		dm_dt[j] = 0;
		dm_dt[k] = 0;
		continue;
       	}

	    double coeff = -gamma / (1 + alpha[index] * alpha[index]);

	    double mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
	    double mh = m[i] * h[i] + m[j] * h[j] + m[k] * h[k];

            //hp=mm.h-mh.m=-mx(mxh)
            double hpi = mm*h[i] - mh*m[i];
            double hpj = mm*h[j] - mh*m[j];
            double hpk = mm*h[k] - mh*m[k];

	    double mth0 = cross_x(m[i], m[j], m[k], hpi, hpj, hpk);
	    double mth1 = cross_y(m[i], m[j], m[k], hpi, hpj, hpk);
	    double mth2 = cross_z(m[i], m[j], m[k], hpi, hpj, hpk);

	    dm_dt[i] = coeff * (mth0 - hpi * alpha[index]);
	    dm_dt[j] = coeff * (mth1 - hpj * alpha[index]);
	    dm_dt[k] = coeff * (mth2 - hpk * alpha[index]);

	    //the above part is standard LLG equation.

	    double coeff_stt = a_J[index] / (1 + alpha[index] * alpha[index]);

	    double mp = m[i] * p[i] + m[j] * p[j] + m[k] * p[k];

	    hpi = mm*p[i] - mp * m[i];
	    hpj = mm*p[j] - mp * m[j];
	    hpk = mm*p[k] - mp * m[k];

	    mth0 = cross_x(m[i], m[j], m[k], hpi, hpj, hpk);
	    mth1 = cross_y(m[i], m[j], m[k], hpi, hpj, hpk);
	    mth2 = cross_z(m[i], m[j], m[k], hpi, hpj, hpk);


	    dm_dt[i] += coeff_stt * ((1 + alpha[index] * beta) * hpi
				     - (beta - alpha[index]) * mth0);
	    dm_dt[j] += coeff_stt * ((1 + alpha[index] * beta) * hpj
				     - (beta - alpha[index]) * mth1);
	    dm_dt[k] += coeff_stt * ((1 + alpha[index] * beta) * hpk
				     - (beta - alpha[index]) * mth2);

	    double c = 6 * sqrt(dm_dt[i] * dm_dt[i] +
				dm_dt[j] * dm_dt[j] +
				dm_dt[k]* dm_dt[k]);

	    dm_dt[i] += c * (1 - mm) * m[i];
	    dm_dt[j] += c * (1 - mm) * m[j];
	    dm_dt[k] += c * (1 - mm) * m[k];

	}

}
