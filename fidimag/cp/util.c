#include "clib.h"

//compute the S \cdot (S_i \times S_j)  
inline double volume(double S[3], double Si[3], double Sj[3]) {
	double tx = S[0] * (-Si[2] * Sj[1] + Si[1] * Sj[2]);
	double ty = S[1] * (Si[2] * Sj[0] - Si[0] * Sj[2]);
	double tz = S[2] * (-Si[1] * Sj[0] + Si[0] * Sj[1]);
	return tx + ty + tz;
}

// C = S_i \dot (S_{i+1} \times S_{j+1}) +  S_i \dot (S_{i-1} \times S_{j-1})
double skyrmion_number(double *spin, double *charge, int nx, int ny, int nz) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j;
	int index, id;
	
	double sum = 0;

	double S[3], S_i[3], S_j[3];

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			index = nyz * i + nz * j;
			S[0] = spin[index];
			S[1] = spin[index + n1];
			S[2] = spin[index + n2];
			
			S_i[0] = S_i[1] = S_i[2] = 0;
			S_j[0] = S_j[1] = S_j[2] = 0;
			if (j > 0) {
				id = index - nz;
				S_j[0] = spin[id];
				S_j[1] = spin[id + n1];
				S_j[2] = spin[id + n2];
			}

			if (i > 0) {
				id = index - nyz;
				S_i[0] = spin[id];
				S_i[1] = spin[id + n1];
				S_i[2] = spin[id + n2];
			}
			
			charge[index]  = volume(S, S_i, S_j);
			
			S_i[0] = S_i[1] = S_i[2] = 0;
			S_j[0] = S_j[1] = S_j[2] = 0;
			if (i < nx - 1 ) {
				id = index + nyz;
				S_i[0] = spin[id];
				S_i[1] = spin[id + n1];
				S_i[2] = spin[id + n2];
			}

			if (j < ny - 1) {
				id = index + nz;
				S_j[0] = spin[id];
				S_j[1] = spin[id + n1];
				S_j[2] = spin[id + n2];
			}


			charge[index]  += volume(S, S_i, S_j);
			charge[index] /= (8*WIDE_PI);

			sum += charge[index];
		}
	}

	return sum;

}

//compute the first derivative respect to x and for the whole mesh
//assume 2d pbc is used 
void compute_px_py_c(double *spin, int nx, int ny, int nz, double *px, double *py){
	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;

	for(int i=0;i<nx; i++){
		for (int j = 0; j < ny; j++) {
			int index = nyz * i + nz * j;


			//x-1
			int id1 = index - nyz;
            if (i==0) {
                id1 += n1;
            }

			//x+1
			int id2 = index + nyz;
            if (i == nx-1){
            	id2 -= n1;
            }

			px[index] = (spin[id2]-spin[id1])/2.0;
			px[index+n1] = (spin[id2 + n1]-spin[id1 + n1])/2.0;
			px[index+n2] = (spin[id2 + n2]-spin[id1 + n2])/2.0;

			//y-1
			id1 = index - nz;
            if (j==0) {
                id1 += nyz;
            }

			//y+1
			id2 = index + nz;
            if (j == ny-1){
                id2 -= nyz;
            }

           	py[index] = (spin[id2]-spin[id1])/2.0;
			py[index+n1] = (spin[id2 + n1]-spin[id1 + n1])/2.0;
			py[index+n2] = (spin[id2 + n2]-spin[id1 + n2])/2.0;
		}
	}

}

// compute the guiding centre, Dynamics of magnetic vortices,     N. Papanicolaou,
// T.N. Tomaras 360, 425-462, (1991)
void compute_guiding_center(double *spin, int nx, int ny, int nz, double *res) {

	int nyz = ny * nz;
	int n1 = nx * nyz, n2 = 2 * n1;
	int i, j;
	int index, id;

	double charge;
	double sum = 0, Rx = 0, Ry = 0;

	double S[3], S_i[3], S_j[3];

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			index = nyz * i + nz * j;
			S[0] = spin[index];
			S[1] = spin[index + n1];
			S[2] = spin[index + n2];

			S_i[0] = S_i[1] = S_i[2] = 0;
			S_j[0] = S_j[1] = S_j[2] = 0;
			if (j > 0) {
				id = index - nz;
				S_j[0] = spin[id];
				S_j[1] = spin[id + n1];
				S_j[2] = spin[id + n2];
			}

			if (i > 0) {
				id = index - nyz;
				S_i[0] = spin[id];
				S_i[1] = spin[id + n1];
				S_i[2] = spin[id + n2];
			}

			charge = volume(S, S_i, S_j);
			sum += charge;
			Rx += i*charge;
			Ry += j*charge;

			S_i[0] = S_i[1] = S_i[2] = 0;
			S_j[0] = S_j[1] = S_j[2] = 0;
			if (i < nx - 1 ) {
				id = index + nyz;
				S_i[0] = spin[id];
				S_i[1] = spin[id + n1];
				S_i[2] = spin[id + n2];
			}

			if (j < ny - 1) {
				id = index + nz;
				S_j[0] = spin[id];
				S_j[1] = spin[id + n1];
				S_j[2] = spin[id + n2];
			}


			charge = volume(S, S_i, S_j);
			sum += charge;
			Rx += i*charge;
			Ry += j*charge;
		}
	}

	res[0]=Rx/sum;
	res[1]=Ry/sum;

}
