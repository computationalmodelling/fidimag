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

void compute_tensors(fft_demag_plan *plan, enum Type_Nij type, double *tensor) {

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

				tensor[id] = NXXdipole(type, x, y, z);
			}
		}
	}
}

fft_demag_plan *create_plan() {

	fft_demag_plan *plan = (fft_demag_plan*) malloc(sizeof(fft_demag_plan));

	return plan;
}

void init_plan(fft_demag_plan *plan, double dx, double dy, double dz, int nx,
		int ny, int nz) {

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
			plan->lenz, plan->tensor_xx, plan->Nxx, FFTW_ESTIMATE);

	plan->m_plan = fftw_plan_dft_r2c_3d(plan->lenx, plan->leny, plan->lenz,
			plan->mx, plan->Mx, FFTW_MEASURE);

	plan->h_plan = fftw_plan_dft_c2r_3d(plan->lenx, plan->leny, plan->lenz,
			plan->Hx, plan->hx, FFTW_MEASURE | FFTW_DESTROY_INPUT);

	printf("hello\n");
	pre_compute(plan);

}

void pre_compute(fft_demag_plan *plan) {
	compute_tensors(plan, Tensor_xx, plan->tensor_xx);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xx, plan->Nyy);

	compute_tensors(plan, Tensor_yy, plan->tensor_yy);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_yy, plan->Nyy);

	compute_tensors(plan, Tensor_zz, plan->tensor_zz);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_zz, plan->Nzz);

	compute_tensors(plan, Tensor_xy, plan->tensor_xy);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xy, plan->Nxy);

	compute_tensors(plan, Tensor_xz, plan->tensor_xz);
	fftw_execute_dft_r2c(plan->tensor_plan, plan->tensor_xz, plan->Nxz);

	compute_tensors(plan, Tensor_yz, plan->tensor_yz);
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

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				id1 = i * nyz + j * nz + k;
				id2 = i * lenyz + j * lenz + k;
				plan->mx[id2] = spin[id1];
				plan->my[id2] = spin[id1 + nxyz];
				plan->mz[id2] = spin[id1 + 2 * nxyz];
			}
		}
	}

	fftw_execute_dft_r2c(plan->m_plan, plan->mx, plan->Mx);
	fftw_execute_dft_r2c(plan->m_plan, plan->my, plan->My);
	fftw_execute_dft_r2c(plan->m_plan, plan->mz, plan->Mz);

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

	for (i = 0; i < plan->total_length; i++) {
		Hx[i] = Nxx[i] * Mx[i] + Nxy[i] * My[i] + Nxz[i] * Mz[i];
		Hy[i] = Nxy[i] * Mx[i] + Nyy[i] * My[i] + Nyz[i] * Mz[i];
		Hz[i] = Nxz[i] * Mx[i] + Nyz[i] * My[i] + Nzz[i] * Mz[i];
	}

	fftw_execute_dft_c2r(plan->m_plan, plan->Hx, plan->hx);
	fftw_execute_dft_c2r(plan->m_plan, plan->Hy, plan->hy);
	fftw_execute_dft_c2r(plan->m_plan, plan->Hz, plan->hz);

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				id1 = i * nyz + j * nz + k;
				id2 = i * lenyz + j * lenz + k;
				field[id1] = plan->hx[id2];
				field[id1 + nxyz] = plan->hy[id2];
				field[id1 + 2 * nxyz] = plan->hz[id2];
			}
		}
	}

}

void print_h(fft_demag_plan *plan){
	int i, j, k, id;

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
				printf("index=%d  hx=%g   hy=%g   hz=%g\n",
						id,
						plan->hx[id],
						plan->hy[id],
						plan->hz[id]);


			}
		}
	}
}


void exact_compute(fft_demag_plan *plan, double *spin, double *field) {
	int i, j, k, i1, i2, i3, index;
	int ii, jj, kk, id1, id2,id3;
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

	for (ii = 0; ii < nx; ii++) {
		for (jj = 0; jj < ny; jj++) {
			for (kk = 0; kk < nz; kk++) {
				id1 = ii * nyz + jj * nz + kk;
				id2 = id1 + nxyz;
				id3 = id2 + nxyz;

				for (i = 0; i < nx; i++) {
					for (j = 0; j < ny; j++) {
						for (k = 0; k < nz; k++) {
							i1 = i * nyz + j * nz + k;
							i2 = i1 + nxyz;
							i3 = i2 + nxyz;
							index = (i - ii + nx - 1) * lenyz + (j - jj + ny
									- 1) * lenz + (k - kk + nz - 1);
							field[id1] += Nxx[index] * spin[i1] + Nxy[index]
									* spin[i2] + Nxz[index] * spin[i3];
							field[id2] += Nxy[index] * spin[i1] + Nyy[index]
									* spin[i2] + Nyz[index] * spin[i3];
							field[id3] += Nxz[index] * spin[i1] + Nyz[index]
									* spin[i2] + Nzz[index] * spin[i3];
						}
					}
				}

			}
		}
	}

}

void finalize_plan(fft_demag_plan *plan) {

	fftw_destroy_plan(plan->m_plan);
	fftw_destroy_plan(plan->h_plan);

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
	printf("at the end\n");
}

