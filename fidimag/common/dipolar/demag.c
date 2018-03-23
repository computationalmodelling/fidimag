#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "dipolar.h"
#include "demagcoef.h"


double Nxxdipole(double x, double y, double z) {
	double x2 = x * x;
	double y2 = y * y;
	double z2 = z * z;
	double R = x2 + y2 + z2;
	if (R == 0)
		return 0.0;
	double r = sqrt(R);
	return -(2 * x2 - y2 - z2) / (R * R * r);
}

double Nxydipole(double x, double y, double z) {
	double R = x * x + y * y + z * z;
	if (R == 0)
		return 0.0;
	double r = sqrt(R);
	return -3 * x * y / (R * R * r);
}

double NXXdipole(enum Type_Nij type, double x, double y, double z) {
	switch (type) {
		case Tensor_xx:
			return Nxxdipole(x, y, z);
		case Tensor_yy:
			return Nxxdipole(y, x, z);
		case Tensor_zz:
			return Nxxdipole(z, y, x);
		case Tensor_xy:
			return Nxydipole(x, y, z);
		case Tensor_xz:
			return Nxydipole(x, z, y);
		case Tensor_yz:
			return Nxydipole(y, z, x);
	}
	return 0;
}

//compute the demag tensors, i.e, H=-N.M
void compute_dipolar_tensors(fft_demag_plan *plan) {

	int i, j, k, id;
	double x, y, z;

	int nx = plan->nx;
	int ny = plan->ny;
	int nz = plan->nz;
	int lenx = plan->lenx;
	int leny = plan->leny;
	int lenz = plan->lenz;
    int lenxy = lenx * leny;

	
	for (k = 0; k < lenz; k++) {
		for (j = 0; j < leny; j++) {
			for (i = 0; i < lenx; i++) {
 				id = k * lenxy + j * lenx + i;
				x = (i - nx + 1) * plan->dx;
				y = (j - ny + 1) * plan->dy;
				z = (k - nz + 1) * plan->dz;

				plan->tensor_xx[id] = NXXdipole(Tensor_xx, x, y, z);
				plan->tensor_yy[id] = NXXdipole(Tensor_yy, x, y, z);
				plan->tensor_zz[id] = NXXdipole(Tensor_zz, x, y, z);
				plan->tensor_xy[id] = NXXdipole(Tensor_xy, x, y, z);
				plan->tensor_xz[id] = NXXdipole(Tensor_xz, x, y, z);
				plan->tensor_yz[id] = NXXdipole(Tensor_yz, x, y, z);

			}
		}
	}
}

//compute the demag tensors, i.e, H=-N.M
void compute_demag_tensors(fft_demag_plan *plan) {

	int i, j, k, id;
	double x, y, z;

	int nx = plan->nx;
	int ny = plan->ny;
	int nz = plan->nz;
	int lenx = plan->lenx;
	int leny = plan->leny;
	int lenz = plan->lenz;
	int lenxy = lenx * leny;

	double dx = plan->dx;
	double dy = plan->dy;
	double dz = plan->dz;

	double length = pow(dx*dy*dz, 1/3.0);
	double asymptotic_radius_sq = pow(26.0*length,2.0);

	for (k = 0; k < lenz; k++) {
		for (j = 0; j < leny; j++) {
			for (i = 0; i < lenx; i++) {
 				id = k * lenxy + j * lenx + i;

				x = (i - nx + 1) * dx;
				y = (j - ny + 1) * dy;
				z = (k - nz + 1) * dz;

				double radius_sq = x*x+y*y+z*z;

				if (radius_sq>asymptotic_radius_sq){
					//printf("%g %g %g %g %g %g\n",x,y,z,dx,dy,dz);
					plan->tensor_xx[id] = DemagNxxAsymptotic(x, y, z, dx, dy, dz);
					plan->tensor_yy[id] = DemagNyyAsymptotic(x, y, z, dx, dy, dz);
					plan->tensor_zz[id] = DemagNzzAsymptotic(x, y, z, dx, dy, dz);
					plan->tensor_xy[id] = DemagNxyAsymptotic(x, y, z, dx, dy, dz);
					plan->tensor_xz[id] = DemagNxzAsymptotic(x, y, z, dx, dy, dz);
					plan->tensor_yz[id] = DemagNyzAsymptotic(x, y, z, dx, dy, dz);
				}else{
					//printf("%g %g %g %g %g %g\n",x,y,z,dx,dy,dz);
					plan->tensor_xx[id] = CalculateSDA00(x, y, z, dx, dy, dz);
					plan->tensor_yy[id] = CalculateSDA11(x, y, z, dx, dy, dz);
					plan->tensor_zz[id] = CalculateSDA22(x, y, z, dx, dy, dz);
					plan->tensor_xy[id] = CalculateSDA01(x, y, z, dx, dy, dz);
					plan->tensor_xz[id] = CalculateSDA02(x, y, z, dx, dy, dz);
					plan->tensor_yz[id] = CalculateSDA12(x, y, z, dx, dy, dz);
				}
			}
		}
	}
}



//used for debug
void print_r(char *str, double *restrict x, int n) {
	int i;
	printf("%s:\n", str);
	for (i = 0; i < n; i++) {
		printf("%f ", x[i]);
	}
	printf("\n");

}

void print_c(char *str, fftw_complex *restrict x, int n) {
	int i;
	printf("%s\n", str);
	for (i = 0; i < n; i++) {
		//printf("%f+%fI  ", x[i]);
	}
	printf("\n");

}

fft_demag_plan *create_plan(void) {

	fft_demag_plan *plan = (fft_demag_plan*) malloc(sizeof(fft_demag_plan));

	return plan;
}

void init_plan(fft_demag_plan *restrict plan, double dx, double dy,
		double dz, int nx, int ny, int nz) {

	//plan->mu_s = mu_s;

	fftw_init_threads();
	fftw_plan_with_nthreads(omp_get_max_threads());

	plan->dx = dx;
	plan->dy = dy;
	plan->dz = dz;

	plan->nx = nx;
	plan->ny = ny;
	plan->nz = nz;

	plan->lenx = 2 * nx - 1;
	plan->leny = 2 * ny - 1;
	plan->lenz = 2 * nz - 1;

	plan->total_length = plan->lenx * plan->leny * plan->lenz;

	int size1 = plan->total_length * sizeof(double);
	int size2 = plan->total_length * sizeof(fftw_complex);

	plan->tensor_xx = (double *) fftw_malloc(size1);
	plan->tensor_yy = (double *) fftw_malloc(size1);
	plan->tensor_zz = (double *) fftw_malloc(size1);
	plan->tensor_xy = (double *) fftw_malloc(size1);
	plan->tensor_xz = (double *) fftw_malloc(size1);
	plan->tensor_yz = (double *) fftw_malloc(size1);

	plan->mx = (double *) fftw_malloc(size1);
	plan->my = (double *) fftw_malloc(size1);
	plan->mz = (double *) fftw_malloc(size1);

	plan->hx = (double *) fftw_malloc(size1);
	plan->hy = (double *) fftw_malloc(size1);
	plan->hz = (double *) fftw_malloc(size1);

	plan->Nxx = (fftw_complex *) fftw_malloc(size2);
	plan->Nyy = (fftw_complex *) fftw_malloc(size2);
	plan->Nzz = (fftw_complex *) fftw_malloc(size2);
	plan->Nxy = (fftw_complex *) fftw_malloc(size2);
	plan->Nxz = (fftw_complex *) fftw_malloc(size2);
	plan->Nyz = (fftw_complex *) fftw_malloc(size2);

	plan->Mx = (fftw_complex *) fftw_malloc(size2);
	plan->My = (fftw_complex *) fftw_malloc(size2);
	plan->Mz = (fftw_complex *) fftw_malloc(size2);
	plan->Hx = (fftw_complex *) fftw_malloc(size2);
	plan->Hy = (fftw_complex *) fftw_malloc(size2);
	plan->Hz = (fftw_complex *) fftw_malloc(size2);

}


void create_fftw_plan(fft_demag_plan *restrict plan) {

	
	plan->tensor_plan = fftw_plan_dft_r2c_3d(plan->lenz, plan->leny,
			plan->lenx, plan->tensor_xx, plan->Nxx,
			FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

	plan->m_plan = fftw_plan_dft_r2c_3d(plan->lenz, plan->leny, plan->lenx,
			plan->mx, plan->Mx, FFTW_MEASURE);

	plan->h_plan = fftw_plan_dft_c2r_3d(plan->lenz, plan->leny, plan->lenx,
			plan->Hx, plan->hx, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	for (int i = 0; i < plan->total_length; i++) {
		plan->Nxx[i] = 0;
		plan->Nyy[i] = 0;
		plan->Nzz[i] = 0;
		plan->Nxy[i] = 0;
		plan->Nxz[i] = 0;
		plan->Nyz[i] = 0;

		plan->mx[i] = 0;
		plan->my[i] = 0;
		plan->mz[i] = 0;

		plan->hx[i] = 0;
		plan->hy[i] = 0;
		plan->hz[i] = 0;
	}


	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xx, plan->Nxx);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_yy, plan->Nyy);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_zz, plan->Nzz);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xy, plan->Nxy);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xz, plan->Nxz);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_yz, plan->Nyz);
	fftw_destroy_plan(plan->tensor_plan);

}





//The computed results doesn't consider the coefficient of \frac{\mu_0}{4 \pi}, the
//reason is in future we can use the following code directly for continuum case
void compute_fields(fft_demag_plan *restrict plan, double *restrict spin, double *restrict mu_s, double *restrict field) {

	int i, j, k, id1, id2;

	int nx = plan->nx;
	int ny = plan->ny;
	int nz = plan->nz;
	int nxy = nx * ny;

	int lenx = plan->lenx;
	int leny = plan->leny;
	int lenxy = lenx * leny;

	for (i = 0; i < plan->total_length; i++) {
		plan->mx[i] = 0;
		plan->my[i] = 0;
		plan->mz[i] = 0;
	}


	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				id1 = k * nxy + j * nx + i;
				id2 = k * lenxy + j * lenx + i;

				plan->mx[id2] = spin[3*id1]*mu_s[id1];
				plan->my[id2] = spin[3*id1+1]*mu_s[id1];
				plan->mz[id2] = spin[3*id1+2]*mu_s[id1];
			}
		}
	}

	//print_r("plan->mx", plan->mx, plan->total_length);

	fftw_execute_dft_r2c(plan->m_plan, plan->mx, plan->Mx);
	fftw_execute_dft_r2c(plan->m_plan, plan->my, plan->My);
	fftw_execute_dft_r2c(plan->m_plan, plan->mz, plan->Mz);

	//print_c("plan->Mx", plan->Mx, plan->total_length);

	fftw_complex *Nxx = plan->Nxx;
	fftw_complex *Nyy = plan->Nyy;
	fftw_complex *Nzz = plan->Nzz;
	fftw_complex *Nxy = plan->Nxy;
	fftw_complex *Nxz = plan->Nxz;
	fftw_complex *Nyz = plan->Nyz;

	fftw_complex *Mx = plan->Mx;
	fftw_complex *My = plan->My;
	fftw_complex *Mz = plan->Mz;
	fftw_complex *Hx = plan->Hx;
	fftw_complex *Hy = plan->Hy;
	fftw_complex *Hz = plan->Hz;

	//print_c("Mx", Mx, plan->total_length);

	for (i = 0; i < plan->total_length; i++) {
		Hx[i] = Nxx[i] * Mx[i] + Nxy[i] * My[i] + Nxz[i] * Mz[i];
		Hy[i] = Nxy[i] * Mx[i] + Nyy[i] * My[i] + Nyz[i] * Mz[i];
		Hz[i] = Nxz[i] * Mx[i] + Nyz[i] * My[i] + Nzz[i] * Mz[i];
	}

	//print_c("Hx", Hx, plan->total_length);

	fftw_execute_dft_c2r(plan->h_plan, plan->Hx, plan->hx);
	fftw_execute_dft_c2r(plan->h_plan, plan->Hy, plan->hy);
	fftw_execute_dft_c2r(plan->h_plan, plan->Hz, plan->hz);
	//print_r("hx", plan->hx, plan->total_length);
	//print_r("hy", plan->hy, plan->total_length);
	//print_r("hz", plan->hz, plan->total_length);

	double scale = -1.0  / plan->total_length;

	for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				id1 = k * nxy + j * nx + i;

				id2 = (k + nz - 1) * lenxy + (j + ny - 1) * lenx + (i + nx - 1);
				field[3*id1] = plan->hx[id2] * scale;
				field[3*id1+1]  = plan->hy[id2] * scale;
				field[3*id1+2]  = plan->hz[id2] * scale;
			}
		}
	}

}

//only used for debug
void exact_compute(fft_demag_plan *restrict plan, double *restrict spin,  double *restrict mu_s, double *restrict field) {
	int i, j, k, index;
	int ip, jp, kp, idf, ids;
	int nx = plan->nx;
	int ny = plan->ny;
	int nz = plan->nz;
        int nxy = nx * ny;

	//int lenx = plan->lenx;
	int lenx = plan->lenx;
	int leny = plan->leny;
        int lenxy = lenx * leny;

	double *Nxx = plan->tensor_xx;
	double *Nyy = plan->tensor_yy;
	double *Nzz = plan->tensor_zz;
	double *Nxy = plan->tensor_xy;
	double *Nxz = plan->tensor_xz;
	double *Nyz = plan->tensor_yz;

	
        for (k = 0; k < nz; k++) {
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				idf = nxy * k + nx * j + i;

				field[3*idf] = 0;
				field[3*idf+1] = 0;
				field[3*idf+2] = 0;

				for (kp = 0; kp < nz; kp++) {
					for (jp = 0; jp < ny; jp++) {
						for (ip = 0; ip < nx; ip++) {
							ids = nxy * kp + nx * jp + ip;
							index = (kp - k + nz - 1) * lenxy + (jp - j + ny
									- 1) * lenx + (ip - i + nx - 1);
							field[3*idf] += (Nxx[index] * spin[3*ids] + Nxy[index]
									* spin[3*ids+1] + Nxz[index] * spin[3*ids+2])*mu_s[ids];
							field[3*idf+1]  += (Nxy[index] * spin[3*ids] + Nyy[index]
									* spin[3*ids+1] + Nyz[index] * spin[3*ids+2])*mu_s[ids];
							field[3*idf+2]  += (Nxz[index] * spin[3*ids] + Nyz[index]
									* spin[3*ids+1] + Nzz[index] * spin[3*ids+2])*mu_s[ids];
						}
					}
				}

				field[3*idf] *= (-1);
				field[3*idf+1]  *= (-1);
				field[3*idf+2]  *= (-1);
			}
		}
	}

}

double compute_demag_energy(fft_demag_plan *restrict plan, double *restrict spin, double *restrict mu_s, double *restrict field) {

	int i,j;

	int nxyz = plan->nx * plan->ny * plan->nz;

	double energy = 0;

	for (i = 0; i < nxyz; i++) {
		j = 3*i;
		energy += mu_s[i]*(spin[j]*field[j]+spin[j+1]*field[j+1]+spin[j+2]*field[j+2]);
	}

	energy = -energy / 2.0;

	return energy;

}

void finalize_plan(fft_demag_plan *restrict plan) {

	fftw_destroy_plan(plan->m_plan);
	fftw_destroy_plan(plan->h_plan);

	fftw_free(plan->tensor_xx);
	fftw_free(plan->tensor_yy);
	fftw_free(plan->tensor_zz);
	fftw_free(plan->tensor_xy);
	fftw_free(plan->tensor_xz);
	fftw_free(plan->tensor_yz);

	fftw_free(plan->Nxx);
	fftw_free(plan->Nyy);
	fftw_free(plan->Nzz);
	fftw_free(plan->Nxy);
	fftw_free(plan->Nxz);
	fftw_free(plan->Nyz);

	fftw_free(plan->Mx);
	fftw_free(plan->My);
	fftw_free(plan->Mz);
	fftw_free(plan->Hx);
	fftw_free(plan->Hy);
	fftw_free(plan->Hz);

	fftw_free(plan->mx);
	fftw_free(plan->my);
	fftw_free(plan->mz);
	fftw_free(plan->hx);
	fftw_free(plan->hy);
	fftw_free(plan->hz);

	fftw_cleanup_threads();

	free(plan);
}

