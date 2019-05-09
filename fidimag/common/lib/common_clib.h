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

inline void normalise(double *m, int *pins, int n){
	int i, j, k;
	double mm;
	for (int id = 0; id < n; id++) {
        i = 3*id;
        j = i + 1;
        k = j + 1;

        if (pins[id]>0) continue;

        mm = sqrt(m[i] * m[i] + m[j] * m[j] + m[k] * m[k]);
        if(mm > 0) {
            mm = 1 / mm;
            m[i] *= mm;
            m[j] *= mm;
            m[k] *= mm;
        }
    }   
}

// ----------------------------------------------------------------------------
// From: llg.c

void llg_rhs(double * dm_dt, double * spin, double * h, double * alpha, int * pins,
             double gamma, int n, int do_precession, double default_c);

void llg_rhs_jtimes(double * jtn, double * m, double * h, double * mp, double * hp,
                    double * alpha, int * pins, double gamma, int n,
                    int do_precession, double default_c);

// ----------------------------------------------------------------------------
// From: stt.c

void compute_stt_field_c(double * spin, double * field, double * jx, double * jy,
                         double * jz, double dx, double dy, double dz, int * ngbs,
                         int n);

void llg_stt_rhs(double * dm_dt, double * m, double * h, double * h_stt,
                 double * alpha, double beta, double u0, double gamma, int n);

void llg_stt_cpp(double * dm_dt, double * m, double * h, double * p, double * alpha,
                 int * pins, double * a_J, double beta, double gamma, int n);

// ----------------------------------------------------------------------------
// From steepest_descent.c

void sd_update_spin (double *spin, double *spin_last, double *magnetisation,
                     double *mxH, double *mxmxH, double *mxmxH_last, double tau,
                     int* pins, int n);

void sd_compute_step (double *spin, double *spin_last, double *magnetisation, double *field, 
                      double *mxH, double *mxmxH, double *mxmxH_last, double tau,
                      int *pins, int n, int counter, double tmin, double tmax);


#endif
