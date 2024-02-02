#include <math.h>

#include <omp.h>

#define WIDE_PI 3.1415926535897932384626433832795L
#define MU0 1.25663706143591728850e-6
#define MU0_INV 795774.71545947669074

inline double cross_x(double a0, double a1, double a2,
                      double b0, double b1, double b2) { return a1 * b2 - a2 * b1; }
inline double cross_y(double a0, double a1, double a2,
                      double b0, double b1, double b2) { return a2 * b0 - a0 * b2; }
inline double cross_z(double a0, double a1, double a2,
                      double b0, double b1, double b2) { return a0 * b1 - a1 * b0; }

void compute_exch_field_micro(double *restrict m, double *restrict field, double *restrict energy, double *restrict Ms_inv,
                              double A, double dx, double dy, double dz, int n, int *ngbs);

void dmi_field(double *restrict m, double *restrict field,
               double *restrict energy, double *restrict Ms_inv,
               double *restrict D, int n_DMIs,
               double *dmi_vector,
               double dx, double dy, double dz, int n, int *ngbs);

void compute_exch_field_rkky_micro(double *m, double *field, double *energy, double *Ms_inv,
                                   double sigma, int nx, double ny, double nz, int z_bottom, int z_top);

void compute_uniaxial_anis(double *restrict m, double *restrict field, double *restrict energy, double *restrict Ms_inv,
                           double *restrict Ku, double *restrict axis, int nx, int ny, int nz);

void compute_uniaxial4_anis(double *restrict m, double *restrict field, double *restrict energy, double *restrict Ms_inv,
                            double *restrict K1, double *restrict K2, double *restrict axis, int nx, int ny, int nz);

double skyrmion_number(double *restrict spin, double *restrict charge,
                       int nx, int ny, int nz, int *restrict ngbs);
