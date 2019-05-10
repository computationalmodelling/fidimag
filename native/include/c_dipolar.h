#pragma once
#include<math.h>
#include<complex>
#include<fftw3.h>
#include "c_vectormath.h"

enum Type_Nij {
	Tensor_xx, Tensor_yy, Tensor_zz, Tensor_xy, Tensor_xz, Tensor_yz
};

//==========================================
//used for demag

typedef struct {
  int nx;
  int ny;
  int nz;
  double dx;
  double dy;
  double dz;
  int lenx;
  int leny;
  int lenz;

  int total_length;

  //TODO: free tensors after obtaining NXX to save memory?
  double *tensor_xx;
  double *tensor_yy;
  double *tensor_zz;
  double *tensor_xy;
  double *tensor_xz;
  double *tensor_yz;

  //TODO: (1)using double, (2)using small size
  std::complex<double> *Nxx; //Nxx, Nxy .. are pure real now.
  std::complex<double> *Nyy;
  std::complex<double> *Nzz;
  std::complex<double> *Nxy;
  std::complex<double> *Nxz;
  std::complex<double> *Nyz;

  std::complex<double> *Mx;
  std::complex<double> *My;
  std::complex<double> *Mz;
  std::complex<double>*Hx;
  std::complex<double> *Hy;
  std::complex<double> *Hz;

  double *mx;
  double *my;
  double *mz;
  double *hx;
  double *hy;
  double *hz;

  //we need three plans
  fftw_plan tensor_plan;
  fftw_plan m_plan;
  fftw_plan h_plan;

} fft_demag_plan;

fft_demag_plan *create_plan(void);
void finalize_plan(fft_demag_plan * plan);
void init_plan(fft_demag_plan *plan, double dx, double dy,
		double dz, int nx, int ny, int nz);
void compute_dipolar_tensors(fft_demag_plan * plan);
void compute_demag_tensors(fft_demag_plan * plan);
void create_fftw_plan(fft_demag_plan * plan);

void compute_demag_tensors_2dpbc(fft_demag_plan * plan, double * tensors, double pbc_2d_error, int sample_repeat_nx, int sample_repeat_ny, double dipolar_radius);
void fill_demag_tensors_c(fft_demag_plan * plan, double * tensors);

void compute_fields(fft_demag_plan * plan, double * spin, double *mu_s, double * field);
void exact_compute(fft_demag_plan * plan, double * spin, double *mu_s, double * field);
double compute_demag_energy(fft_demag_plan * plan, double * spin, double * mu_s, double * field, double * energy);



//=========================================================
//=========================================================
//used for sode
typedef struct {
	int nxyz;

	double dt;
	double T;
	double gamma;
	double *mu_s;
	double coeff;
	double Q;

	double theta;
	double theta1;
	double theta2;

	double *dm1;
	double *dm2;
	double *eta;

} ode_solver;

void init_solver(ode_solver *s, double k_B, double theta, int nxyz, double dt, double gamma);
ode_solver *create_ode_plan(void);
void finalize_ode_plan(ode_solver *plan);
void run_step1(ode_solver *s, double *m, double *h, double *m_pred, double *T,
		double *alpha, double *mu_s_inv, int *pins);
void run_step2(ode_solver *s, double *m_pred, double *h, double *m, double *T,
		double *alpha, double *mu_s_inv, int *pins);
