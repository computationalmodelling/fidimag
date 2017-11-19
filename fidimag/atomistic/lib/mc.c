#include "clib.h"
#include "fidimag_random.h"

inline double dot(double a[3], double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// helper function to compute the bulk DMI for a cuboid mesh with neighbour i
inline double dmi_energy_site(double a[3], double b[3], int i) {
  double res = 0;
  switch (i) {
  case 0:
    res = -(a[1] * b[2] - a[2] * b[1]);
    break; //-x
  case 1:
    res = (a[1] * b[2] - a[2] * b[1]);
    break; //+x
  case 2:
    res = -(a[2] * b[0] - a[0] * b[2]);
    break; //-y
  case 3:
    res = (a[2] * b[0] - a[0] * b[2]);
    break; //+y
  case 4:
    res = -(a[0] * b[1] - a[1] * b[0]);
    break; //-z
  case 5:
    res = (a[0] * b[1] - a[1] * b[0]);
    break; //+z
  default:
    break;
  }
  return res;
}

// note that this DMI only works for a 2d hexagnoal mesh.
inline double dmi_energy_hexagnoal_site(double a[3], double b[3], int i) {
  double res = 0;
  switch (i) {
  case 0:
    res = (a[1] * b[2] - a[2] * b[1]);
    break; //+x
  case 1:
    res = -(a[1] * b[2] - a[2] * b[1]);
    break; //-x
  case 2:
    res = 0.5 * (a[1] * b[2] - a[2] * b[1]) +
          0.86602540378443864676 * (a[2] * b[0] - a[0] * b[2]);
    break; // top right, r={1/2, sqrt(3)/2, 0}
  case 3:
    res = -0.5 * (a[1] * b[2] - a[2] * b[1]) -
          0.86602540378443864676 * (a[2] * b[0] - a[0] * b[2]);
    break; // bottom left, r={-1/2, -sqrt(3)/2, 0}
  case 4:
    res = -0.5 * (a[1] * b[2] - a[2] * b[1]) +
          0.86602540378443864676 * (a[2] * b[0] - a[0] * b[2]);
    break; // top left, r={-1/2, sqrt(3)/2, 0}
  case 5:
    res = 0.5 * (a[1] * b[2] - a[2] * b[1]) -
          0.86602540378443864676 * (a[2] * b[0] - a[0] * b[2]);
    break; // bottom right, r={1/2, -sqrt(3)/2, 0}
  default:
    break;
  }
  return res;
}

// helper function to compute the cubic energy
inline double cubic_energy_site(double *m, double Kc) {

  double mx2 = m[0] * m[0];
  double my2 = m[1] * m[1];
  double mz2 = m[2] * m[2];

  return -Kc * (mx2 * mx2 + my2 * my2 + mz2 * mz2);
}

double compute_deltaE_anisotropy(double *spin, double *new_spin, double *h,
                                 double Kc, int i, int new_i) {

  double energy1 = -dot(&spin[3 * i], &h[3 * i]); // zeeman energy
  double energy2 = -dot(&new_spin[3 * new_i], &h[3 * i]);

  energy1 += cubic_energy_site(&spin[3 * i], Kc); // cubic anisotropy energy
  energy2 += cubic_energy_site(&new_spin[3 * new_i], Kc);

  return energy2 - energy1;
}

/*
 * n is the total spin number
 */
double compute_deltaE_exchange_DMI_hexagnoal(double *spin, double *new_spin,
                                             int *ngbs, double J, double D,
                                             int i, int new_i) {

  int id_nn = 6 * i;
  double energy1 = 0, energy2 = 0;

  for (int j = 0; j < 6; j++) {
    int k = ngbs[id_nn + j];
    if (k >= 0) {
      energy1 -= J * dot(&spin[3 * i], &spin[3 * k]); // exchange energy
      energy1 += D * dmi_energy_hexagnoal_site(&spin[3 * i], &spin[3 * k],
                                               j); // DMI energy

      energy2 -= J * dot(&new_spin[3 * new_i], &spin[3 * k]);
      energy2 += D * dmi_energy_hexagnoal_site(&new_spin[3 * new_i],
                                               &spin[3 * k], j); // DMI energy
    }
  }

  return energy2 - energy1;
}

double compute_deltaE_exchange_DMI(double *spin, double *new_spin, int *ngbs,
                                   int *nngbs, double J, double J1, double D,
                                   double D1, int i, int new_i) {

  int id_nn = 6 * i;
  double energy1 = 0, energy2 = 0;

  for (int j = 0; j < 6; j++) {
    int k = ngbs[id_nn + j];
    if (k >= 0) {
      energy1 -= J * dot(&spin[3 * i], &spin[3 * k]); // exchange energy
      energy1 += D * dmi_energy_site(&spin[3 * i], &spin[3 * k], j); // DMI

      energy2 -= J * dot(&new_spin[3 * new_i], &spin[3 * k]);
      energy2 += D * dmi_energy_site(&new_spin[3 * new_i], &spin[3 * k], j);
    }
    k = nngbs[id_nn + j];
    if (k >= 0) {
      energy1 -= J1 * dot(&spin[3 * i], &spin[3 * k]); // exchange energy
      energy1 += D1 * dmi_energy_site(&spin[3 * i], &spin[3 * k], j); // DMI

      energy2 -= J1 * dot(&new_spin[3 * new_i], &spin[3 * k]);
      energy2 += D1 * dmi_energy_site(&new_spin[3 * new_i], &spin[3 * k], j);
    }
  }

  return energy2 - energy1;
}

void run_step_mc(mt19937_state *state, double *spin, double *new_spin,
                 int *ngbs, int *nngbs, double J, double J1, double D,
                 double D1, double *h, double Kc, int n, double T,
                 int hexagnoal_mesh) {

  double delta_E, r;
  int update = 0;

  uniform_random_sphere(state, new_spin, n);

  for (int new_i = 0; new_i < n; new_i++) {

    int i = rand_int_n(state, n);
    int j = 3 * i;
    int new_j = 3 * new_i;

    delta_E =
        compute_deltaE_anisotropy(&spin[0], &new_spin[0], &h[0], Kc, i, new_i);

    if (hexagnoal_mesh) {
      delta_E += compute_deltaE_exchange_DMI_hexagnoal(
          &spin[0], &new_spin[0], &ngbs[0], J, D, i, new_i);
    } else {
      delta_E += compute_deltaE_exchange_DMI(&spin[0], &new_spin[0], &ngbs[0],
                                             &nngbs[0], J, J1, D, D1, i, new_i);
    }

    update = 0;

    if (delta_E < 0) {
      update = 1;
    } else {
      r = random_double_half_open(state);
      if (r < exp(-delta_E / T))
        update = 1;
    }

    if (update) {
      spin[j] = new_spin[new_j];
      spin[j + 1] = new_spin[new_j + 1];
      spin[j + 2] = new_spin[new_j + 2];
    }
  }
}
