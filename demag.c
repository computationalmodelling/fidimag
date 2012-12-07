#include <math.h>
#include "clib.h"

enum Type_Nij {
	xx, yy, zz, xy, xz, yz
};

double Nxxdipole(double x, double y, double z) {
	double x2 = x * x;
	double y2 = y * y;
	double z2 = z * z;
	double R = x2 + y2 + z2;
	double r = sqrt(R);
	return -(2 * x2 - y2 - z2) / (R * R * r);
}

double Nxydipole(double x, double y, double z) {
	double R = x * x + y * y + z * z;
	double r = sqrt(R);
	return -3 * x * y / (R * R * r);
}

double NXXdipole(enum Type_Nij type, double x, double y, double z) {
	switch (type) {
	case xx:
		return Nxxdipole(x, y, z);
	case yy:
		return Nxxdipole(y, x, z);
	case zz:
		return Nxxdipole(z, y, x);
	case xy:
		return Nxydipole(x, y, z);
	case xz:
		return Nxydipole(x, z, y);
	case yz:
		return Nxydipole(y, z, x);
	}
	return 0;
}

void compute_tensors(double *Nxx, double *Nyy, double *Nzz, double *Nxy,
		double *Nxz, double *Nyz, double dx, double dy, double dz, int nx,
		int ny, int nz) {

	int i, j, k, id;
	double x, y, z;

	int lenx = 2 * nx - 1;
	int leny = 2 * ny - 1;
	int lenz = 2 * nz - 1;
	int lenyz = leny * lenz;

	for (i = 0; i < lenx; i++) {
		for (j = 0; j < leny; j++) {
			for (k = 0; k < lenz; k++) {
				id = i * lenyz + j * lenz + k;
				x = (i - nx + 1) * dx;
				y = (j - ny + 1) * dy;
				z = (k - nz + 1) * dz;

				Nxx[id] = NXXdipole(xx, x, y, z);
				Nyy[id] = NXXdipole(yy, x, y, z);
				Nzz[id] = NXXdipole(zz, x, y, z);
				Nxy[id] = NXXdipole(xy, x, y, z);
				Nxx[id] = NXXdipole(xz, x, y, z);
				Nyz[id] = NXXdipole(yz, x, y, z);
			}
		}
	}
}


