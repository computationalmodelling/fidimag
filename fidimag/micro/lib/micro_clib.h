#include<math.h>

//#include<omp.h>

#define WIDE_PI 3.1415926535897932384626433832795L
#define MU0 1.25663706143591728850e-6
#define MU0_INV 795774.71545947669074

inline double cross_x(double a0, double a1, double a2, double b0, double b1, double b2) { return a1*b2 - a2*b1; }
inline double cross_y(double a0, double a1, double a2, double b0, double b1, double b2) { return a2*b0 - a0*b2; }
inline double cross_z(double a0, double a1, double a2, double b0, double b1, double b2) { return a0*b1 - a1*b0; }

void compute_exch_field_micro(double *m, double *field, double *energy, double *Ms_inv,
                         double A, double dx, double dy, double dz, int nx, int ny, int nz, int xperiodic, int yperiodic, int zperiodic);

void dmi_field_bulk(double *m, double *field, double *energy, double *Ms_inv,
                    double *D, double dx, double dy, double dz,
                    int nx, int ny, int nz, int xperiodic, int yperiodic, int zperiodic);

void dmi_field_interfacial(double *m, double *field, double *energy, double *Ms_inv,
                           double *D, double dx, double dy, double dz,
                           int nx, int ny, int nz, int xperiodic, int yperiodic);


void compute_uniaxial_anis(double *m, double *field, double *energy, double *Ms_inv, 
	double *Ku, double *axis, int nx, int ny, int nz);
