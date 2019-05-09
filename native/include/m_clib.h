#include<math.h>
#include<omp.h>

#define WIDE_PI 3.1415926535897932384626433832795L
#define MU0 1.25663706143591728850e-6
#define MU0_INV 795774.71545947669074

double cross_x(double a0, double a1, double a2, double b0, double b1, double b2);
double cross_y(double a0, double a1, double a2, double b0, double b1, double b2);
double cross_z(double a0, double a1, double a2, double b0, double b1, double b2);

void compute_exch_field_micro(double * m, double * field, double * energy, double * Ms_inv,
                         double A, double dx, double dy, double dz, int n, int *ngbs);

void dmi_field(double * m, double * field,
               double * energy, double * Ms_inv,
               double * D, int n_DMIs,
               double *dmi_vector,
               double dx, double dy, double dz, int n, int *ngbs);

void compute_exch_field_rkky_micro(double *m, double *field, double *energy, double *Ms_inv,
                         double sigma, int nx, double ny, double nz, int z_bottom, int z_top);

void compute_uniaxial_anis(double * m, double * field, double * energy, double * Ms_inv,
	double * Ku, double * axis, int nx, int ny, int nz);

double skyrmion_number(double * spin, double * charge,
                       int nx, int ny, int nz, int * ngbs);
