#ifndef __CLIB__
#define __CLIB__

#include <math.h>
//#include<omp.h>
#define WIDE_PI 3.1415926535897932384626433832795L

// ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------
// From: llg.c

void llg_rhs(double *restrict dm_dt, double *restrict spin, double *restrict h, double *restrict alpha, int *restrict pins,
             double gamma, int n, int do_precession, double default_c);

void llg_rhs_jtimes(double *restrict jtn, double *restrict m, double *restrict h, double *restrict mp, double *restrict hp,
                    double *restrict alpha, int *restrict pins, double gamma, int n,
                    int do_precession, double default_c);

// ----------------------------------------------------------------------------
// From: stt.c

void compute_stt_field_c(double *restrict spin, double *restrict field, double *restrict jx, double *restrict jy,
                         double *restrict jz, double dx, double dy, double dz, int *restrict ngbs,
                         int n);

void llg_stt_rhs(double *restrict dm_dt, double *restrict m, double *restrict h, double *restrict h_stt,
                 double *restrict alpha, double beta, double u0, double gamma, int n);

void llg_stt_cpp(double *restrict dm_dt, double *restrict m, double *restrict h, double *restrict p, double *restrict alpha,
                 int *restrict pins, double *restrict a_J, double beta, double gamma, int n);


#endif
