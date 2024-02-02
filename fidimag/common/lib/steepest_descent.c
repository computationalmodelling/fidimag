#include "common_clib.h"
#include "math.h"

void sd_update_spin(double *spin, double *spin_last, double *magnetisation,
                    double *mxH, double *mxmxH, double *mxmxH_last, double tau,
                    int *pins, int n) {

// Update the field -------------------------------------------------------
#pragma omp parallel for
    for (int i = 0; i < n; i++) {

        if (magnetisation[i] > 0.0) {

            int spin_idx;
            double mxH_sq;
            double factor_plus;
            double factor_minus;
            double new_spin[3];

            spin_idx = 3 * i;

            // Copy spin from previous (counter) step
            for (int j = 0; j < 3; j++) {
                spin_last[spin_idx + j] = spin[spin_idx + j];
            }

            mxH_sq = 0;
            for (int j = 0; j < 3; j++) {
                mxH_sq += mxH[spin_idx + j] * mxH[spin_idx + j];
            }

            factor_plus = 4 + tau * tau * mxH_sq;
            factor_minus = 4 - tau * tau * mxH_sq;

            for (int j = 0; j < 3; j++) {
                new_spin[j] = factor_minus * spin[spin_idx + j] - 4 * tau * mxmxH[spin_idx + j];
                new_spin[j] = new_spin[j] / factor_plus;

                spin[spin_idx + j] = new_spin[j];
            }
        } // close if Ms or mu_s > 0
    }     // close for
    normalise(spin, pins, n);
}

void sd_compute_step(double *spin, double *spin_last, double *magnetisation, double *field,
                     double *mxH, double *mxmxH, double *mxmxH_last, double tau,
                     int *pins, int n, int counter, double tmin, double tmax) {

    double res;
    double num = 0;
    double den = 0;
    double sign;
#pragma omp parallel for reduction(+ : num, den)
    for (int i = 0; i < n; i++) {
        if (magnetisation[i] > 0.0) {
            int spin_idx;
            double ds[3], dy[3];

            spin_idx = 3 * i;

            // Copy mxmxH from previous (counter) step
            for (int j = 0; j < 3; j++) {
                mxmxH_last[spin_idx + j] = mxmxH[spin_idx + j];
            }

            // Compute the torques

            mxH[spin_idx] = cross_x(spin[spin_idx], spin[spin_idx + 1], spin[spin_idx + 2],
                                    field[spin_idx],
                                    field[spin_idx + 1],
                                    field[spin_idx + 2]);

            mxH[spin_idx + 1] = cross_y(spin[spin_idx], spin[spin_idx + 1], spin[spin_idx + 2],
                                        field[spin_idx],
                                        field[spin_idx + 1],
                                        field[spin_idx + 2]);

            mxH[spin_idx + 2] = cross_z(spin[spin_idx], spin[spin_idx + 1], spin[spin_idx + 2],
                                        field[spin_idx],
                                        field[spin_idx + 1],
                                        field[spin_idx + 2]);

            mxmxH[spin_idx] = cross_x(spin[spin_idx], spin[spin_idx + 1], spin[spin_idx + 2],
                                      mxH[spin_idx], mxH[spin_idx + 1], mxH[spin_idx + 2]);
            mxmxH[spin_idx + 1] = cross_y(spin[spin_idx], spin[spin_idx + 1], spin[spin_idx + 2],
                                          mxH[spin_idx], mxH[spin_idx + 1], mxH[spin_idx + 2]);
            mxmxH[spin_idx + 2] = cross_z(spin[spin_idx], spin[spin_idx + 1], spin[spin_idx + 2],
                                          mxH[spin_idx], mxH[spin_idx + 1], mxH[spin_idx + 2]);

            // Compute terms for the calculation of the time step tau
            for (int j = 0; j < 3; j++) {
                ds[j] = spin[spin_idx + j] - spin_last[spin_idx + j];
                dy[j] = mxmxH[spin_idx + j] - mxmxH_last[spin_idx + j];
            }

            if (counter % 2 == 0) {
                for (int j = 0; j < 3; j++) {
                    num += ds[j] * ds[j];
                    den += ds[j] * dy[j];
                }
            } else {
                for (int j = 0; j < 3; j++) {
                    num += ds[j] * dy[j];
                    den += dy[j] * dy[j];
                }
            }

        } // Close if Ms or mu_s > 0
    }     // close for

    // Criteria for the evaluation of tau is in line 96 of:
    // https://github.com/MicroMagnum/MicroMagnum/blob/minimizer/src/magnum/micromagnetics/micro_magnetics_solver.py
    if (den == 0.0) {
        res = tmax;
    } else {
        res = num / den;
    }

    // Micromagnum criteria:
    sign = (res > 0) ? 1 : ((res < 0) ? -1 : 0);
    tau = fmax(fmin(fabs(res), tmax), tmin) * sign;

    // In MuMax3 they only define a specific value:
    // https://github.com/mumax/3/blob/master/engine/minimizer.go
    // tau = res;
}
