#ifndef __CLIB__
#define __CLIB__

#include <complex.h>
#include <fftw3.h>
#include <math.h>
//#include<omp.h>

#include "fidimag_random.h"

#define WIDE_PI 3.1415926535897932384626433832795L

/* 3 components for the cross product calculations */
inline double cross_x(double a0, double a1, double a2, double b0, double b1,
                      double b2) {
  return a1 * b2 - a2 * b1;
}
inline double cross_y(double a0, double a1, double a2, double b0, double b1,
                      double b2) {
  return a2 * b0 - a0 * b2;
}
inline double cross_z(double a0, double a1, double a2, double b0, double b1,
                      double b2) {
  return a0 * b1 - a1 * b0;
}

void compute_exch_field(double *spin, double *field, double *energy, double Jx,
                        double Jy, double Jz, int *ngbs, int n, int n_ngbs);

void compute_exch_field_spatial(double *spin, double *field, double *energy,
                                double *J, int *ngbs, int n, int n_ngbs);

double compute_exch_energy(double *spin, double Jx, double Jy, double Jz,
                           int nx, int ny, int nz, int xperiodic,
                           int yperiodic);

void compute_anis(double *spin, double *field, double *energy, double *Ku,
                  double *axis, int n);
void compute_anis_cubic(double *spin, double *field, double *energy,
	          double *Kc, int n); 

void dmi_field_bulk(double *spin, double *field, double *energy, double *D,
                    int *ngbs, int n, int n_ngbs);

void dmi_field_interfacial_atomistic(double *spin, double *field,
                                     double *energy, double D, int *ngbs, int n,
                                     int n_ngbs, int n_ngbs_dmi, double *DMI_vec);

void demag_full(double *spin, double *field, double *energy, double *coords,
                double *mu_s, double *mu_s_scale, int n);

double dmi_energy(double *spin, double D, int nx, int ny, int nz, int xperiodic,
                  int yperiodic);

void llg_s_rhs(double *dm_dt, double *spin, double *h, double *alpha,
               double *chi, double gamma, int n);

void normalise(double *m, int *pins, int n);

double skyrmion_number(double *spin, double *charge, int nx, int ny, int nz,
                       int *ngbs, int n_ngbs);

double skyrmion_number_BergLuscher(double *spin, double *charge, int nx, int ny,
                                   int nz, int *ngbs, int n_ngbs);

void compute_guiding_center(double *spin, int nx, int ny, int nz, int nx_start,
                            int nx_stop, int ny_start, int ny_stop,
                            double *res);

void compute_px_py_c(double *spin, int nx, int ny, int nz, double *px,
                     double *py);

//======================================================================

void llg_rhs_dw_c(double *m, double *h, double *dm, double *T, double *alpha,
                  double *mu_s_inv, int *pins, double *eta, int n, double gamma,
                  double dt);

//======================================================================

void run_step_mc(mt19937_state *state, double *spin, double *new_spin,
                 int *ngbs, int *nngbs, double J, double J1, double D,
                 double D1, double *h, double Kc, int n, double T,
                 int hexagnoal_mesh);

#endif
