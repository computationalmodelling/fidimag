#include "clib.h"

/* Classical mean-field (CMF) relaxation of a lattice of classical spins, after Hovorka and Sluckin (arXiv:2007.12777).
 *
 * In the mean-field approximation each spin is a spin in the effective field of its neighbours' averages, and the thermal average of a
 * classical spin of unit length in the field B (tesla) with moment mu_s is
 *
 *     <S> = L(x) B / |B|,    x = mu_s |B| / (k_B T),    L(x) = coth(x) - 1/x         (the Langevin function)
 *
 * The effective field is that of the energy terms of the simulation evaluated at the average spins, so the kernels here do not know the
 * Hamiltonian: they take the field that has been computed and give the average it implies. The self-consistent solution is found by
 * iterating S <- S + alpha (target(S) - S), or with Anderson mixing. Spins that are pinned, or that have no moment, are not changed.
 */

/* L(x) / x, which is what turns the field into the moment: <S> = (L(x) / x) mu_s B / k_B T. It tends to 1/3 for a vanishing field, so the
 * direction of the field is never needed and a field of zero is no special case. The series is used for small x, where coth(x) - 1/x
 * would lose digits. */
static inline double langevin_over_x(double x) {
    if (x < 0.1) {
        double x2 = x * x;
        return 1.0 / 3 - x2 / 45 + 2 * x2 * x2 / 945 - x2 * x2 * x2 / 4725;
    }
    if (x > 30) return (1.0 - 1.0 / x) / x;
    return (1.0 / tanh(x) - 1.0 / x) / x;
}

/* resid = target(S) - S for the field `field`, and the largest absolute value of a component of it. T is the temperature in kelvin, k_B the
 * Boltzmann constant of the same units as mu_s B (joules); at T = 0 the target is the unit vector along the field, and a site where the
 * field is zero is left as it is. */
double cmf_residual(double *restrict spin, double *restrict field, double *restrict mu_s, int *restrict pins,
                    double *restrict resid, int n, double T, double k_B) {
    double max_resid = 0;

#pragma omp parallel for reduction(max : max_resid)
    for (int i = 0; i < n; i++) {
        int j = 3 * i;
        if (pins[i] > 0 || mu_s[i] <= 0) {
            resid[j] = resid[j + 1] = resid[j + 2] = 0;
            continue;
        }
        double B = sqrt(field[j] * field[j] + field[j + 1] * field[j + 1] + field[j + 2] * field[j + 2]);
        double f;
        if (T > 0) {
            double beta_mu = mu_s[i] / (k_B * T);
            f = langevin_over_x(beta_mu * B) * beta_mu;
        } else {
            f = B > 0 ? 1.0 / B : 0.0;
        }
        for (int c = 0; c < 3; c++) {
            double r = (B > 0 || T > 0) ? f * field[j + c] - spin[j + c] : 0.0;
            resid[j + c] = r;
            if (fabs(r) > max_resid) max_resid = fabs(r);
        }
    }
    return max_resid;
}

/* The linear mixing step: S += alpha resid (the residual is zero where a spin does not change) */
void cmf_mix(double *restrict spin, double *restrict resid, double alpha, int n3) {
#pragma omp parallel for
    for (int i = 0; i < n3; i++) spin[i] += alpha * resid[i];
}

#define CMF_MAX_HISTORY 32

/* Anderson mixing. X holds nh earlier spin arrays and R their residuals (row k at X + k * 3n and R + k * 3n, the last one the current). The
 * combination of the residuals of smallest norm with weights that add up to 1, c = A^-1 1 / (1 A^-1 1) for the matrix A of their scalar
 * products, gives the new spins sum_k c_k (X_k + alpha R_k), which are then cut to a length of at most 1. A is made a little positive
 * definite (as the residuals become parallel when it converges), and if it cannot be solved the plain step from the last row is taken.
 * Returns 0, or 1 for the plain step, or -1 if there are too many rows. */
int cmf_anderson(double *restrict spin, double *restrict X, double *restrict R, int nh, double alpha, double *restrict mu_s,
                 int *restrict pins, int n) {
    if (nh < 1 || nh > CMF_MAX_HISTORY) return -1;
    long n3 = 3L * n;
    double A[CMF_MAX_HISTORY][CMF_MAX_HISTORY], b[CMF_MAX_HISTORY], c[CMF_MAX_HISTORY];
    int plain = 0;

    for (int k = 0; k < nh; k++) {
        for (int l = 0; l <= k; l++) {
            double s = 0;
            const double *rk = R + k * n3, *rl = R + l * n3;
#pragma omp parallel for reduction(+ : s)
            for (long i = 0; i < n3; i++) s += rk[i] * rl[i];
            A[k][l] = A[l][k] = s;
        }
    }
    double trace = 0;
    for (int k = 0; k < nh; k++) trace += A[k][k];
    if (nh == 1 || !(trace > 0)) {
        plain = 1;
    } else {
        for (int k = 0; k < nh; k++) {
            A[k][k] += 1e-10 * trace / nh;
            b[k] = 1.0;
        }
        /* Gaussian elimination with partial pivoting */
        for (int col = 0; col < nh && !plain; col++) {
            int piv = col;
            for (int r = col + 1; r < nh; r++)
                if (fabs(A[r][col]) > fabs(A[piv][col])) piv = r;
            if (!(fabs(A[piv][col]) > 1e-300)) { plain = 1; break; }
            if (piv != col) {
                for (int q = 0; q < nh; q++) { double t = A[col][q]; A[col][q] = A[piv][q]; A[piv][q] = t; }
                double t = b[col]; b[col] = b[piv]; b[piv] = t;
            }
            for (int r = col + 1; r < nh; r++) {
                double m = A[r][col] / A[col][col];
                for (int q = col; q < nh; q++) A[r][q] -= m * A[col][q];
                b[r] -= m * b[col];
            }
        }
        if (!plain) {
            double total = 0;
            for (int r = nh - 1; r >= 0; r--) {
                double s = b[r];
                for (int q = r + 1; q < nh; q++) s -= A[r][q] * c[q];
                c[r] = s / A[r][r];
                total += c[r];
            }
            if (!(fabs(total) > 1e-300) || !isfinite(total)) plain = 1;
            else for (int k = 0; k < nh; k++) c[k] /= total;
            for (int k = 0; k < nh; k++) if (!isfinite(c[k])) plain = 1;
        }
    }
    if (plain) {
        for (int k = 0; k < nh; k++) c[k] = 0;
        c[nh - 1] = 1;
    }

#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        if (pins[i] > 0 || mu_s[i] <= 0) continue;
        double v[3];
        for (int a = 0; a < 3; a++) {
            long j = 3L * i + a;
            double s = 0;
            for (int k = 0; k < nh; k++) s += c[k] * (X[k * n3 + j] + alpha * R[k * n3 + j]);
            v[a] = s;
        }
        double len = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        double scale = len > 1.0 ? 1.0 / len : 1.0;
        for (int a = 0; a < 3; a++) spin[3L * i + a] = v[a] * scale;
    }
    return plain;
}
