#include <math.h>
#include <stdlib.h>
#include "dipolar.h"
#include "demagcoef.h"
#include "tensor_2dpbc.h"

double DemagTensorNormal(enum Type_Nij comp,double x,double y,double z,double a,double b,double c)
{
	switch(comp){
		case Tensor_xx:return CalculateSDA00(x,y,z,a,b,c);
		case Tensor_yy:return CalculateSDA00(y,x,z,b,a,c);
		case Tensor_zz:return CalculateSDA00(z,y,x,c,b,a);
		case Tensor_xy:return CalculateSDA01(x,y,z,a,b,c);
		case Tensor_xz:return CalculateSDA01(x,z,y,a,c,b);
		case Tensor_yz:return CalculateSDA01(y,z,x,b,c,a);
	}

  return -1;
}


double DemagTensorAsymptotic(enum Type_Nij comp,double x,double y,double z,double a,double b,double c)
{
	switch(comp){
		case Tensor_xx:return DemagNxxAsymptotic(x,y,z,a,b,c);
		case Tensor_yy:return DemagNxxAsymptotic(y,x,z,b,a,c);
		case Tensor_zz:return DemagNxxAsymptotic(z,y,x,c,b,a);
		case Tensor_xy:return DemagNxyAsymptotic(x,y,z,a,b,c);
		case Tensor_xz:return DemagNxyAsymptotic(x,z,y,a,c,b);
		case Tensor_yz:return DemagNxyAsymptotic(y,z,x,b,c,a);
	}
 
  return -1;

}

double  DemagTensorDipolar(enum Type_Nij comp,double x,double y,double z){
   
     switch(comp){
	case Tensor_xx:return DemagNxxDipolar(x,y,z);
        case Tensor_yy:return DemagNxxDipolar(y,x,z);
	case Tensor_zz:return DemagNxxDipolar(z,y,x);
	case Tensor_xy:return DemagNxyDipolar(x,y,z);
	case Tensor_xz:return DemagNxyDipolar(x,z,y);
	case Tensor_yz:return DemagNxyDipolar(y,z,x);
    }

  return -1;

}


double DemagTensorInfinite(enum Type_Nij comp, double x,double y,double z,double X0,double Y0)
{
	switch(comp){
		case Tensor_xx:return Nxxinf(x,y,z,X0,Y0);
		case Tensor_yy:return Nyyinf(x,y,z,X0,Y0);
		case Tensor_zz:return Nzzinf(x,y,z,X0,Y0);
		case Tensor_xy:return Nxyinf(x,y,z,X0,Y0);
		case Tensor_xz:return Nxzinf(x,y,z,X0,Y0);
		case Tensor_yz:return Nyzinf(x,y,z,X0,Y0);
	}

  return -1;
}


double compute_single_tensor(enum Type_Nij comp, int g, double x, double y, double z,
        double a, double b, double c, int xdim, int ydim, 
		double asymptotic_radius_sq, double dipolar_radius_sq){

    if ((comp == Tensor_xy || comp == Tensor_xz || comp == Tensor_yz) && x * y == 0) return 0.0;

    double Tx = xdim*a, Ty = ydim*b, cof1 = 1 / (4 * WIDE_PI * a * b * c), cof2 = a * b * c / (4 * WIDE_PI);
    double* tmpx = (double *)malloc(sizeof(double)*(2 * g + 2));
    double* tmpy = (double *)malloc(sizeof(double)*(2 * g + 2));

    double tpx, tpy, radius_sq;
    for (int i = -g; i <= g; i++) {
        for (int j = -g; j <= g; j++) {
            tpx = x + i*Tx;
            tpy = y + j*Ty;
            radius_sq = tpx * tpx + tpy * tpy + z*z;
            if (radius_sq < asymptotic_radius_sq) {
                tmpy[j + g] = DemagTensorNormal(comp, tpx, tpy, z, a, b, c) * cof1;
            } else if (radius_sq < dipolar_radius_sq) {
                tmpy[j + g] = DemagTensorAsymptotic(comp, tpx, tpy, z, a, b, c);
            } else {
                //printf("%f\n", radius_sq);
                tmpy[j + g] = DemagTensorDipolar(comp, tpx, tpy, z) * cof2;
            }
        }
        tmpx[i + g] = AccurateSum(2 * g + 1, tmpy);
    }

    double X0 = (g + 0.5) * Tx;
    double Y0 = (g + 0.5) * Ty;

    tmpx[2 * g + 1] = DemagTensorInfinite(comp, x, y, z, X0, Y0) * cof2 / (Tx * Ty);

    double result = AccurateSum(2 * g + 2, tmpx);

    free(tmpx);
    free(tmpy);

    return result;
}

double compute_single_tensor_finitely(enum Type_Nij comp, int gx, int gy, double x, double y, double z,
        double a, double b, double c, int xdim, int ydim, 
		double asymptotic_radius_sq, double dipolar_radius_sq){

    if ((comp == Tensor_xy || comp == Tensor_xz || comp == Tensor_yz) && x * y == 0) return 0.0;

    double Tx = xdim*a, Ty = ydim*b, cof1 = 1 / (4 * WIDE_PI * a * b * c), cof2 = a * b * c / (4 * WIDE_PI);
    double* tmpx = (double *)malloc(sizeof(double)*(2 * gx + 2));
    double* tmpy = (double *)malloc(sizeof(double)*(2 * gy + 2));
    double tpx, tpy, radius_sq;
    for (int i = -gx; i <= gx; i++) {
        for (int j = -gy; j <= gy; j++) {
            tpx = x + i*Tx;
            tpy = y + j*Ty;
            radius_sq = tpx * tpx + tpy * tpy + z*z;
            if (radius_sq < asymptotic_radius_sq) {
                tmpy[j + gy] = DemagTensorNormal(comp, tpx, tpy, z, a, b, c) * cof1;
            } else if (radius_sq < dipolar_radius_sq) {
                tmpy[j + gy] = DemagTensorAsymptotic(comp, tpx, tpy, z, a, b, c);
            } else {
                // printf("%f\n", radius_sq);
                tmpy[j + gy] = DemagTensorDipolar(comp, tpx, tpy, z) * cof2;
            }
        }
        tmpx[i + gx] = AccurateSum(2 * gy + 1, tmpy);
    }


    double result = AccurateSum(2 * gx + 1, tmpx);

    free(tmpx);
    free(tmpy);

    return result;
}


int FindG(enum Type_Nij comp, double v, double Tx, double Ty, double pbc_2d_error){

    double tmp;

    switch (comp) {
        case Tensor_xy:
        case Tensor_xz:
        case Tensor_yz:
        case Tensor_xx:
            tmp = v / (4 * WIDE_PI * (Tx * Tx) * sqrt(Tx * Tx + Ty * Ty) * pbc_2d_error);
            return (int) pow(tmp, 1 / 3.0) + 1;
        case Tensor_yy:
            tmp = v / (4 * WIDE_PI * (Ty * Ty) * sqrt(Tx * Tx + Ty * Ty) * pbc_2d_error);
            return (int) pow(tmp, 1 / 3.0) + 1;
        case Tensor_zz:
            tmp = v * sqrt(Tx * Tx + Ty * Ty) / (4 * WIDE_PI * (Tx * Ty * Tx * Ty) * pbc_2d_error);
            return (int) pow(tmp, 1 / 3.0) + 1;
    }
  return 0;
}



//compute the demag tensors, i.e, H=-N.M
void compute_demag_tensors_2dpbc(fft_demag_plan *plan, double pbc_2d_error, int sample_repeat_nx, int sample_repeat_ny, double dipolar_radius){
   
      int gxx = sample_repeat_nx;
      int gyy = sample_repeat_ny;
      int gzz = 0;

	int i, j, k, id;
	double x, y, z;

	int nx = plan->nx;
	int ny = plan->ny;
	int nz = plan->nz;
	int lenx = plan->lenx;
	int leny = plan->leny;
	int lenz = plan->lenz;
	int len_xy = lenx * leny;

	double dx = plan->dx;
	double dy = plan->dy;
	double dz = plan->dz;

        double maxedge = dx;
    	if (dy > maxedge) maxedge = dy;
    	if (dz > maxedge) maxedge = dz;
    	dx /= maxedge;
    	dy /= maxedge;
   	dz /= maxedge;

	double length = pow(dx*dy*dz, 1/3.0);
	double asymptotic_radius_sq = pow(26.0*length,2.0);
	double dipolar_radius_sq = dipolar_radius*dipolar_radius;

       if (gxx >= 0 && gyy >= 0) {

	  for (k = 0; k < lenz; k++) {
		for (j = 0; j < leny; j++) {
			for (i = 0; i < lenx; i++) {
 				id = k * len_xy + j * lenx + i;

				x = (i - nx + 1) * dx;
				y = (j - ny + 1) * dy;
				z = (k - nz + 1) * dz;

				plan->tensor_xx[id] = compute_single_tensor_finitely(Tensor_xx, gxx, gyy, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_yy[id] = compute_single_tensor_finitely(Tensor_yy, gxx, gyy, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_zz[id] = compute_single_tensor_finitely(Tensor_zz, gxx, gyy, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_xy[id] = compute_single_tensor_finitely(Tensor_xy, gxx, gyy, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_xz[id] = compute_single_tensor_finitely(Tensor_xz, gxx, gyy, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_yz[id] = compute_single_tensor_finitely(Tensor_yz, gxx, gyy, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
			}
		}
	}
     }else{


        gxx = FindG(Tensor_xx, dx * dy*dz, nx*dx, ny * dy, pbc_2d_error);
        gyy = FindG(Tensor_yy, dx * dy*dz, nx*dx, ny * dy, pbc_2d_error);
        gzz = FindG(Tensor_zz, dx * dy*dz, nx*dx, ny * dy, pbc_2d_error);

	   for (k = 0; k < lenz; k++) {
		for (j = 0; j < leny; j++) {
			for (i = 0; i < lenx; i++) {
 				id = k * len_xy + j * lenx + i;

				x = (i - nx + 1) * dx;
				y = (j - ny + 1) * dy;
				z = (k - nz + 1) * dz;

				plan->tensor_xx[id] = compute_single_tensor(Tensor_xx, gxx, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_yy[id] = compute_single_tensor(Tensor_yy, gyy, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_zz[id] = compute_single_tensor(Tensor_zz, gzz, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_xy[id] = compute_single_tensor(Tensor_xy, gxx, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_xz[id] = compute_single_tensor(Tensor_xz, gxx, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
				plan->tensor_yz[id] = compute_single_tensor(Tensor_yz, gxx, x, y, z, dx, dy, dz, nx, ny, asymptotic_radius_sq, dipolar_radius_sq);
			}
		}
	  }

   }
	
}





