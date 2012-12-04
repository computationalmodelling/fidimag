
void compute_uniform_exch(double *spin,double *field,double J,double dx,double dy,double dz,int nx,int ny, int nz);
void compute_anis(double *spin,double *field,double Kx,double Ky,double Kz,int nx,int ny, int nz);
void llg_rhs(double *dm_dt,double *spin,double *h,double gamma,double alpha, double mu_s, int nxyz);
