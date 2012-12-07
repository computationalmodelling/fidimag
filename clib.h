void compute_uniform_exch(double *spin, double *field, double J, double dx,
		double dy, double dz, int nx, int ny, int nz);

void compute_anis(double *spin, double *field, double Dx, double Dy, double Dz,
		int nxyz);

void llg_rhs(double *dm_dt, double *spin, double *h, double gamma,
		double alpha, double mu_s, int nxyz, double c);

void compute_tensors(double *Nxx, double *Nyy, double *Nzz, double *Nxy,
		double *Nxz, double *Nyz, double dx, double dy, double dz, int nx,
		int ny, int nz);
