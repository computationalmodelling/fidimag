#include<math.h>
#include<complex.h>
#include<fftw3.h>
//#include<omp.h>

#define WIDE_PI 3.1415926535897932384626433832795L

/* 3 components for the cross product calculations */
inline double cross_x(double a0, double a1, double a2,
                      double b0, double b1, double b2) { return a1 * b2 - a2 * b1; }
inline double cross_y(double a0, double a1, double a2,
                      double b0, double b1, double b2) { return a2 * b0 - a0 * b2; }
inline double cross_z(double a0, double a1, double a2,
                      double b0, double b1, double b2) { return a0 * b1 - a1 * b0; }

void compute_exch_field(double *spin, double *field, double *energy,
						double Jx, double Jy, double Jz,
                        int *ngbs, int n);
double compute_exch_energy(double *spin, double Jx, double Jy, double Jz,
                           int nx, int ny, int nz,
                           int xperiodic, int yperiodic);

void compute_anis(double *spin, double *field, double *energy,
	              double *Ku, double *axis, int n);


void dmi_field_bulk(double *spin, double *field, double *energy,
                    double D, int *ngbs, int n);

void dmi_field_interfacial_atomistic(double *spin, double *field,
                                     double *energy, double D, int *ngbs,
                                     int n, int nneighbours,
                                     double *DMI_vec);

void demag_full(double *spin, double *field, double *energy,
                double *coords, double *mu_s, int n)

double dmi_energy(double *spin, double D, int nx, int ny, int nz,
                  int xperiodic, int yperiodic);

void llg_rhs(double * dm_dt, double * spin, double * h, double *alpha,
		int *pins, double gamma, int n, int do_procession, double default_c);

void llg_rhs_jtimes(double *jtn, double *m, double *h,
                    double *mp, double *hp, double *alpha, int *pins,
                    double gamma, int n, int do_procession, double default_c);

void llg_s_rhs(double * dm_dt, double * spin, double * h, double *alpha,
             double *chi, double gamma, int n);


void compute_stt_field_c(double *spin, double *field, double *jx, double *jy,
		double dx, double dy, int *ngbs, int n);

void llg_stt_rhs(double *dm_dt, double *m, double *h,
                 double *h_stt, double *alpha,
                 double beta, double u0, double gamma, int n);



void normalise(double *m, int n);

double skyrmion_number(double *spin, double *charge,
                       int nx, int ny, int nz, int *ngbs);

void compute_guiding_center(double *spin, int nx, int ny, int nz, double *res);

void compute_px_py_c(double *spin, int nx, int ny, int nz,
                     double *px, double *py);


//=========================================================
//=========================================================
//used for sode
typedef struct {
	int n;

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

void init_solver(ode_solver *s, double k_B, double theta,
                 int n, double dt, double gamma);

ode_solver *create_ode_plan(void);

void finalize_ode_plan(ode_solver *plan);

void run_step1(ode_solver *s, double *m, double *h, double *m_pred, double *T,
		double *alpha, double *mu_s_inv, int *pins);

void run_step2(ode_solver *s, double *m_pred, double *h, double *m, double *T,
		double *alpha, double *mu_s_inv, int *pins);
