#include<math.h>

//#include<omp.h>

#define WIDE_PI 3.1415926535897932384626433832795L
#define MU0 1.25663706143591728850e-6

inline double cross_x(double a0, double a1, double a2, double b0, double b1, double b2) { return a1*b2 - a2*b1; }
inline double cross_y(double a0, double a1, double a2, double b0, double b1, double b2) { return a2*b0 - a0*b2; }
inline double cross_z(double a0, double a1, double a2, double b0, double b1, double b2) { return a0*b1 - a1*b0; }

void compute_laplace_m(double *m, double *field, double *Ms, double dx, double dy, double dz,
		int nx, int ny, int nz);

void compute_relaxation_field_c(double *m, double *field, double *Ms, double chi_inv, int n);

void compute_perp_field_c(double *m, double *field, double *field_p, int n);

void llg_rhs_baryakhtar(double *dm_dt, double *m, double *h, double *delta_h,
		double *alpha, double beta, int *pins,
		double gamma, int nxyz, int do_precession);


void llg_rhs_baryakhtar_reduced(double *dm_dt, double *m, double *hp, double *delta_hp,
                                double *alpha, double beta, int *pins,
                                double gamma, int nxyz, int do_precession, double default_c);
