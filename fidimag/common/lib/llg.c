#include "common_clib.h"

/* The right hand side of the LLG equation for the CVOde solver.  This can be
 * used both for the micromagnetic and atomistic codes since m or S are unitless
 * and the prefactors keep the same structure.
 *
 * The LLG equation has the structure:
 * ( * is for dot or scalar product)
 *
 *      dm        -gamma
 *     ---- =    --------  ( m X H_eff  + a * m X ( m x H_eff ) )
 *      dt             2
 *              ( 1 + a  )
 *
 * where a is the Gilbert damping constant, gamma is the gyromagnetic
 * ratio ( gamma = 1.76e11 for a free electron; for the micromagnetic
 * case, we use gamma_0 = mu_0 * gamma = 2.21e5, for a free electron),
 * m is the magnetisation vector and H_eff is the effective field.
 *
 * In our calculation, we usually compute:

 *      m x (m x H_eff) = ( m * H_eff) m - (m * m) H_eff

 * Then, if we define the perpendicular component of the effective
 * field as:

 *      H_perp = H_eff - ( m * H_eff) m

 * we can write

 *      dm        -gamma
 *     ---- =    --------  ( m X H_perp - a * H_perp )
 *      dt             2
 *              ( 1 + a  )

 * since m X m = 0 and, theoretically, m * m = 1 (at zero temperature m has
 * fixed length). However, for the second term to the right hand side,
 * proportional to a, it is better to keep the m X ( m x H_eff ) term with (m *
 * m) for stability purposes, thus we use
 *

 *      H_perp = (m * m) H_eff - ( m * H_eff) m

 * for both terms.
 *
 * Additionally, to preserve the magnetisation length, we need to correct the
 * dm / dt term for every time step, adding the following term to the right
 * hand size expression of the LLG:

 *      dm              dm               2
 *     ----    --->    ---- + c * ( 1 - m  ) m
 *      dt              dt

 *  with
                    _________
 *                  /        2
 *      c = 6 *    / (  dm  )
 *                /  ( ---- )
 *              \/   (  dt  )

 * The correction must be introduced since numerically, m can change length
 * during a step of the integration (which would not occur if the integration
 * step is infinitely small), deviating from the real answer. If we just
 * rescaled m, we would have to recompute the effective field (m changes) and
 * also the solution would present jumps due to the rescaling. With the term
 * specified above, m stays close to 1, giving a more continuous solution, and
 * also the term stays as zero if m is exactly 1 (notice that if m increases,
 * the correction is negative and it decreases m length; similarly if m
 * decreases). The prefactor c must be specified because the correction term
 * must be sufficiently strong to affect the solution. Accordingly, we can
 * think of dm/dt as a kind of velocity that is proportional to the change of
 * rate of m, hence using its magnitude, the correction is stronger for large
 * deviations and weak for small deviations.
 *
 * What c does, and how to choose it
 * ---------------------------------
 *
 * There is no physics in this term. It is a numerical stabilisation of the
 * invariant |m| = 1, of the kind Baumgarte introduced for constrained
 * mechanics. Writing rho = |m|, and using that the transverse field above is
 * exactly perpendicular to m, so that the LLG part contributes nothing to
 * d|m|/dt,
 *
 *      d rho / dt = c * rho * ( 1 - rho^2 )
 *
 * which is the normal form of the supercritical pitchfork bifurcation:
 * rho = 0 is unstable and rho = 1 is stable. Linearising about rho = 1 with
 * rho = 1 + eps gives
 *
 *      d eps / dt = -2 c eps
 *
 * so c is nothing more than the decay rate of length errors, with time
 * constant 1 / (2c). That is the whole content of the term, and it is what
 * sets the two competing limits: too small and errors are not removed within
 * a step, too large and the equation is stiffened by an eigenvalue -2c.
 *
 * The factor 6 in c = 6 |dm/dt| therefore means that length errors relax
 * twelve times faster than the magnetisation rotates, which is why it is
 * written against |dm/dt| rather than as an absolute rate: it is scale free.
 * The micromagnetic driver instead defaults to a fixed 1e11, i.e. a time
 * constant of 5 ps, against a precession period of about 35 ps at gamma Ms.
 *
 * Measured on standard problem 4, against the OOMMF reference in
 * examples/micromagnetic/std4. The deviation quoted is the largest difference,
 * over all sample times and the three components, between the spatially
 * averaged magnetisation and the reference: a worst case along the
 * trajectory, not an average, and an error of the averaged magnetisation
 * rather than a per cell one. The drift is the largest ||m_i| - 1| over every
 * cell and sample time.
 *
 *  - At rtol = atol = 1e-10 the choice makes no difference to the answer at
 *    all. Sweeping c over 0, 1e9 ... 1e13 and 6|dm/dt| gives a deviation of
 *    3.20e-05 in every case, including with the term switched off. Only the
 *    length drift responds, falling as 1/c from 5.8e-09 to 5.2e-11, and at
 *    1e13 the step count rises by 15%, which is the stiffening above.
 *  - At rtol = atol = 1e-6, where it earns its keep, the term nearly halves
 *    the error: 6.36e-05 with c = 0 against 3.58e-05 at 1e11.
 *  - For the explicit ARKODE methods it does nothing useful at any tolerance,
 *    and at 1e-8 it makes the drift worse, 3.0e-10 with c = 0 against 1.8e-08
 *    with 6|dm/dt|, since the term is itself a perturbation. Those integrators
 *    can project instead; see the normalise option of ErkSolver.
 *
 * So the value is not an accuracy knob to be tuned. Anything from about 1e11
 * to 1e12, or the 6|dm/dt| rule, behaves the same; it only has to be large
 * enough to remove drift within a step and small enough not to stiffen the
 * equation.
 *
 * The same term, and the same pitchfork analysis, are in A. E. Botha,
 * Stabilisation of the Landau-Lifshitz-Gilbert equation for numerical solution
 * via standard methods, Sci. Rep. 15, 15775 (2025),
 * doi:10.1038/s41598-025-99966-x, which argues the coefficient should be unity
 * once the equation is in a suitable dimensionless form. That paper attributes
 * the |dm/dt| rule used here to Sec. 2.4 of D. Chernyshenko, Computational
 * Methods in Micromagnetics, PhD thesis, University of Southampton (2016).
 */

void llg_rhs(double *restrict dm_dt, double *restrict m, double *restrict h, double *restrict alpha, int *restrict pins,
             double gamma, int n, int do_precession, double default_c) {

    int i, j, k;

    double coeff, mm, mh, c;
    double hpi, hpj, hpk;

#pragma omp parallel for private(i, j, k, coeff, mm, mh, c, hpi, hpj, hpk)
    for (int id = 0; id < n; id++) {
        // Indexes for the 3 components of the spin (magnetic moment)
        // at the i-th lattice (mesh) site  --> x, y, z
        i = 3 * id;
        j = i + 1;
        k = i + 2;

        // Pinned spins do not follow the dynamical equation
        if (pins[id] > 0) {
            dm_dt[i] = 0;
            dm_dt[j] = 0;
            dm_dt[k] = 0;
            continue;
        }

        coeff = -gamma / (1.0 + alpha[id] * alpha[id]);

        // Dot products
        mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
        mh = m[i] * h[i] + m[j] * h[j] + m[k] * h[k];

        // Usually, m is normalised, i.e., mm=1;
        // so hp = mm.h - mh.m = -m x (m x h)
        // We set here the perpendicular componenet of the field
        // but using the (m * m) product
        hpi = mm * h[i] - mh * m[i];
        hpj = mm * h[j] - mh * m[j];
        hpk = mm * h[k] - mh * m[k];

        // Keep mm rather than assuming it is one. The point is that
        //
        //     (m.m) h - (m.h) m
        //
        // is perpendicular to m exactly, whatever the length of m, while
        //
        //     h - (m.h) m
        //
        // leaves a component along m of (m.h)(1 - m.m). That vanishes only
        // when |m| is exactly one, and the correction term below holds |m|
        // near one rather than at one, so the residual is not zero and feeds
        // a spurious change of the spin length back into the dynamics.
        //
        // A note here used to say that dropping mm cost two decimals of
        // accuracy against OOMMF on standard problem 4 (2016). That no longer
        // reproduces. Measured against the OOMMF reference in
        // examples/micromagnetic/std4, over its full 1 ns and 1000 sample
        // points, the deviation of the mean magnetisation is
        //
        //     with mm     max 3.5294e-05   rms 1.3139e-05
        //     without mm  max 3.5515e-05   rms 1.3221e-05
        //
        // a difference of 0.6%, and the two runs differ from each other by
        // 2.2e-07 at most. Whatever was seen in 2016 has been fixed
        // elsewhere since. The mm form is kept because it is the correct
        // projection and costs one multiplication
        double mth0 = 0, mth1 = 0, mth2 = 0;

        // The first term: m x H_eff = m x H_perp
        if (do_precession) {
            mth0 = cross_x(m[i], m[j], m[k], hpi, hpj, hpk);
            mth1 = cross_y(m[i], m[j], m[k], hpi, hpj, hpk);
            mth2 = cross_z(m[i], m[j], m[k], hpi, hpj, hpk);
        }

        // The RHS term of the LLG equation
        dm_dt[i] = coeff * (mth0 - hpi * alpha[id]);
        dm_dt[j] = coeff * (mth1 - hpj * alpha[id]);
        dm_dt[k] = coeff * (mth2 - hpk * alpha[id]);

        // In future, we will try the new method to integrate the LLG equation,
        // A mixed mid-point Runge-Kutta like scheme for the integration of
        // Landau-Lifshitz equation Journal of Applied Physics 115, 17D101
        // (2014) if possible, we can combine it with adaptive step size, don't
        // know how to do but it's worth a try.

        if (default_c < 0) {
            c = 6 * sqrt(dm_dt[i] * dm_dt[i] +
                         dm_dt[j] * dm_dt[j] +
                         dm_dt[k] * dm_dt[k]);
        } else {
            c = default_c;
        }
        // printf("%0.15g   %0.15g\n", c, default_c);

        // Correct the RHS term to keep m normalised
        dm_dt[i] += c * (1 - mm) * m[i];
        dm_dt[j] += c * (1 - mm) * m[j];
        dm_dt[k] += c * (1 - mm) * m[k];
    }
}

void llg_rhs_jtimes(double *restrict jtn, double *restrict m, double *restrict h, double *restrict mp, double *restrict hp, double *restrict alpha, int *restrict pins,
                    double gamma, int n, int do_precession, double default_c) {

    // #pragma omp parallel for private(i,j,k,coeff,mm, mh, c, hpi,hpj,hpk)
    for (int id = 0; id < n; id++) {
        int i = 3 * id;
        int j = i + 1;
        int k = i + 2;

        if (pins[i] > 0) {
            continue;
        }

        double coeff = -gamma / (1.0 + alpha[i] * alpha[i]);

        if (do_precession) {
            jtn[i] = coeff * (cross_x(mp[i], mp[j], mp[k], h[i], h[j], h[k]) + cross_x(m[i], m[j], m[k], hp[i], hp[j], hp[k]));
            jtn[j] = coeff * (cross_y(mp[i], mp[j], mp[k], h[i], h[j], h[k]) + cross_y(m[i], m[j], m[k], hp[i], hp[j], hp[k]));
            jtn[k] = coeff * (cross_z(mp[i], mp[j], mp[k], h[i], h[j], h[k]) + cross_z(m[i], m[j], m[k], hp[i], hp[j], hp[k]));
        } else {
            jtn[i] = 0;
            jtn[j] = 0;
            jtn[k] = 0;
        }

        double mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
        double mh = m[i] * h[i] + m[j] * h[j] + m[k] * h[k];
        double mhp = m[i] * hp[i] + m[j] * hp[j] + m[k] * hp[k];
        double mph = mp[i] * h[i] + mp[j] * h[j] + mp[k] * h[k];
        double mmp = m[i] * mp[i] + m[j] * mp[j] + m[k] * mp[k];

        jtn[i] += alpha[i] * coeff * ((mph + mhp) * m[i] + mh * mp[i] - 2 * mmp * h[i] - mm * hp[i]);
        jtn[j] += alpha[i] * coeff * ((mph + mhp) * m[j] + mh * mp[j] - 2 * mmp * h[j] - mm * hp[j]);
        jtn[k] += alpha[i] * coeff * ((mph + mhp) * m[k] + mh * mp[k] - 2 * mmp * h[k] - mm * hp[k]);

        if (default_c > 0) {
            jtn[i] += default_c * ((1 - mm) * mp[i] - 2 * mmp * m[i]);
            jtn[j] += default_c * ((1 - mm) * mp[j] - 2 * mmp * m[j]);
            jtn[k] += default_c * ((1 - mm) * mp[k] - 2 * mmp * m[k]);
        }
    }
}
