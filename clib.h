#include <complex.h>
#include <fftw3.h>

void compute_uniform_exch(double *spin, double *field, double J, double dx,
		double dy, double dz, int nx, int ny, int nz);

void compute_anis(double *spin, double *field, double Dx, double Dy, double Dz,
		int nxyz);

void llg_rhs(double *dm_dt, double *spin, double *h, double gamma,
		double alpha, double mu_s, int nxyz, double c);

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
	double *tensor;
	double *tensor_xx;
	double *tensor_yy;
	double *tensor_zz;
	double *tensor_xy;
	double *tensor_xz;
	double *tensor_yz;

	fftw_complex *Nxx;
	fftw_complex *Nyy;
	fftw_complex *Nzz;
	fftw_complex *Nxy;
	fftw_complex *Nxz;
	fftw_complex *Nyz;

	fftw_complex *Mx;
	fftw_complex *My;
	fftw_complex *Mz;
	fftw_complex *Hx;
	fftw_complex *Hy;
	fftw_complex *Hz;
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

fft_demag_plan *create_plan();
void finalize_plan(fft_demag_plan *plan);
void init_plan(fft_demag_plan *plan, double dx, double dy, double dz, int nx,
		int ny, int nz);
void compute_fields(fft_demag_plan *plan, double *spin, double *field);
void exact_compute(fft_demag_plan *plan, double *spin, double *field);
void print_h(fft_demag_plan *plan);

