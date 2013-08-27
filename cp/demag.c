#include <math.h>
#include <stdlib.h>
#include "clib.h"

enum Type_Nij {
	Tensor_xx, Tensor_yy, Tensor_zz, Tensor_xy, Tensor_xz, Tensor_yz
};

double Nxxdipole(double x, double y, double z) {
	double x2 = x * x;
	double y2 = y * y;
	double z2 = z * z;
	double R = x2 + y2 + z2;
	if (R == 0)
		return 0.0;
	double r = sqrt(R);
	return (2 * x2 - y2 - z2) / (R * R * r);
}

double Nxydipole(double x, double y, double z) {
	double R = x * x + y * y + z * z;
	if (R == 0)
		return 0.0;
	double r = sqrt(R);
	return 3 * x * y / (R * R * r);
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
void compute_all_tensors(fft_demag_plan *plan) {

	int i, j, k, id;
	double x, y, z;

	int nx = plan->nx;
	int ny = plan->ny;
	int nz = plan->nz;
	int lenx = plan->lenx;
	int leny = plan->leny;
	int lenz = plan->lenz;
	int lenyz = leny * lenz;

	for (i = 0; i < lenx; i++) {
		for (j = 0; j < leny; j++) {
			for (k = 0; k < lenz; k++) {
				id = i * lenyz + j * lenz + k;
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

//used for debug
void print_r(char *str, double *x, int n) {
	int i;
	printf("%s:\n", str);
	for (i = 0; i < n; i++) {
		printf("%f ", x[i]);
	}
	printf("\n");

}

void print_c(char *str, fftw_complex *x, int n) {
	int i;
	printf("%s\n", str);
	for (i = 0; i < n; i++) {
		printf("%f+%fI  ", x[i]);
	}
	printf("\n");

}

fft_demag_plan *create_plan() {

	fft_demag_plan *plan = (fft_demag_plan*) malloc(sizeof(fft_demag_plan));

	return plan;
}

void init_plan(fft_demag_plan *plan, double mu_s, double dx, double dy,
		double dz, int nx, int ny, int nz) {

	plan->mu_s = mu_s;

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

	plan->tensor_plan = fftw_plan_dft_r2c_3d(plan->lenx, plan->leny,
			plan->lenz, plan->tensor_xx, plan->Nxx,
			FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

	plan->m_plan = fftw_plan_dft_r2c_3d(plan->lenx, plan->leny, plan->lenz,
			plan->mx, plan->Mx, FFTW_MEASURE);

	plan->h_plan = fftw_plan_dft_c2r_3d(plan->lenx, plan->leny, plan->lenz,
			plan->Hx, plan->hx, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	int i;
	for (i = 0; i < plan->total_length; i++) {
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

	compute_all_tensors(plan);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xx, plan->Nxx);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_yy, plan->Nyy);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_zz, plan->Nzz);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xy, plan->Nxy);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xz, plan->Nxz);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_yz, plan->Nyz);
	fftw_destroy_plan(plan->tensor_plan);

}

void compute_fields(fft_demag_plan *plan, double *spin, double *field) {

	int i, j, k, id1, id2;

	int nx = plan->nx;
	int ny = plan->ny;
	int nz = plan->nz;
	int nyz = ny * nz;
	int nxyz = nx * nyz;

	int leny = plan->leny;
	int lenz = plan->lenz;
	int lenyz = leny * lenz;

	for (i = 0; i < plan->total_length; i++) {
		plan->mx[i] = 0;
		plan->my[i] = 0;
		plan->mz[i] = 0;
	}

	double *spin_x = &spin[0];
	double *spin_y = &spin[nxyz];
	double *spin_z = &spin[2 * nxyz];
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				id1 = i * nyz + j * nz + k;
				id2 = i * lenyz + j * lenz + k;
				plan->mx[id2] = spin_x[id1];
				plan->my[id2] = spin_y[id1];
				plan->mz[id2] = spin_z[id1];
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

	double scale = plan->mu_s  / plan->total_length;
	double *field_x = &field[0];
	double *field_y = &field[nxyz];
	double *field_z = &field[2 * nxyz];
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				id1 = i * nyz + j * nz + k;
				id2 = (i + nx - 1) * lenyz + (j + ny - 1) * lenz + (k + nz - 1);
				field_x[id1] = plan->hx[id2] * scale;
				field_y[id1] = plan->hy[id2] * scale;
				field_z[id1] = plan->hz[id2] * scale;
			}
		}
	}

}

//only used for debug
void exact_compute(fft_demag_plan *plan, double *spin, double *field) {
	int i, j, k, index;
	int ip, jp, kp, idf, ids;
	int nx = plan->nx;
	int ny = plan->ny;
	int nz = plan->nz;
	int nyz = ny * nz;
	int nxyz = nx * nyz;

	int lenx = plan->lenx;
	int leny = plan->leny;
	int lenz = plan->lenz;
	int lenyz = leny * lenz;

	double *Nxx = plan->tensor_xx;
	double *Nyy = plan->tensor_yy;
	double *Nzz = plan->tensor_zz;
	double *Nxy = plan->tensor_xy;
	double *Nxz = plan->tensor_xz;
	double *Nyz = plan->tensor_yz;

	double *f_x = &field[0];
	double *f_y = &field[nxyz];
	double *f_z = &field[2 * nxyz];

	double *s_x = &spin[0];
	double *s_y = &spin[nxyz];
	double *s_z = &spin[2 * nxyz];
	double scale = plan->mu_s;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				idf = i * nyz + j * nz + k;

				f_x[idf] = 0;
				f_y[idf] = 0;
				f_z[idf] = 0;

				for (ip = 0; ip < nx; ip++) {
					for (jp = 0; jp < ny; jp++) {
						for (kp = 0; kp < nz; kp++) {
							ids = ip * nyz + jp * nz + kp;
							index = (i - ip + nx - 1) * lenyz + (j - jp + ny
									- 1) * lenz + (k - kp + nz - 1);
							f_x[idf] += Nxx[index] * s_x[ids] + Nxy[index]
									* s_y[ids] + Nxz[index] * s_z[ids];
							f_y[idf] += Nxy[index] * s_x[ids] + Nyy[index]
									* s_y[ids] + Nyz[index] * s_z[ids];
							f_z[idf] += Nxz[index] * s_x[ids] + Nyz[index]
									* s_y[ids] + Nzz[index] * s_z[ids];
						}
					}
				}

				f_x[idf] *= scale;
				f_y[idf] *= scale;
				f_z[idf] *= scale;

			}
		}
	}

}

double compute_demag_energy(fft_demag_plan *plan, double *spin, double *field) {

	int i;

	int nxyz = plan->nx * plan->ny * plan->nz;

	double energy = 0, mu_s = plan->mu_s;

	for (i = 0; i < 3 * nxyz; i++) {

		energy += mu_s * spin[i] * field[i];
	}

	energy = -energy / 2.0;

	return energy;

}

void finalize_plan(fft_demag_plan *plan) {

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

	free(plan);
}

