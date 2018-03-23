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
 * deviations and weak for small deviations.  The factor 6 is added ad-hoc,
* which seems to work well when computing the solutions, but its specification
* stills requires a more strict proof.  It is worth mentioning that the norm of
* dm/dt changes the time scaling by a factor proportional to 1/t, therefore in
* the future we could try to estimate its
 * influence with more mathematical/numerical rigour and analyse an optimal
 * value for the prefactor (6). 
 */

void llg_rhs(double *restrict dm_dt, double *restrict m, double *restrict h, double *restrict alpha, int *restrict pins,
        double gamma, int n, int do_precession, double default_c) {

    int i, j, k;

    double coeff, mm, mh, c;
    double hpi, hpj, hpk;

    #pragma omp parallel for private(i,j,k,coeff,mm, mh, c, hpi,hpj,hpk)
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

        // IMPORTANT: do not ignore mm !!
        // What we've found is that if we igonre mm, i.e. using
        //    hpi = h[i] - mh * m[i];
        //    hpj = h[j] - mh * m[j];
        //    hpk = h[k] - mh * m[k];
        // the micromagnetic standard problem 4 failed to converge (?)
        //
        // NOTE (Fri 08 Jul 2016 13:58): In fact, the problem converges but with 2 less
        // decimals of accuracy, compared with the OOMMF calculation
        double mth0 = 0, mth1 = 0, mth2 = 0;

        // The first term: m x H_eff = m x H_perp
        if (do_precession){
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

        if (default_c < 0){
            c = 6 * sqrt(dm_dt[i] * dm_dt[i] +
                         dm_dt[j] * dm_dt[j] + 
                         dm_dt[k] * dm_dt[k]
                         );
        } else {
            c = default_c;
        }
        //printf("%0.15g   %0.15g\n", c, default_c);

        // Correct the RHS term to keep m normalised
        dm_dt[i] += c * (1 - mm) * m[i];
        dm_dt[j] += c * (1 - mm) * m[j];
        dm_dt[k] += c * (1 - mm) * m[k];

    }

}


void llg_rhs_jtimes(double *restrict jtn, double *restrict m, double *restrict h, double *restrict mp, double *restrict hp, double *restrict alpha, int *restrict pins,
        double gamma, int n, int do_precession, double default_c) {

    //#pragma omp parallel for private(i,j,k,coeff,mm, mh, c, hpi,hpj,hpk)
    for (int id = 0; id < n; id++) {
        int i = 3 * id;
        int j = i + 1;
        int k = i + 2;


        if (pins[i]>0){
             continue;
        }

        double coeff = -gamma/(1.0+alpha[i]*alpha[i]);

        if (do_precession){
            jtn[i] = coeff*(cross_x(mp[i],mp[j],mp[k],h[i],h[j],h[k])+cross_x(m[i],m[j],m[k],hp[i],hp[j],hp[k]));
            jtn[j] = coeff*(cross_y(mp[i],mp[j],mp[k],h[i],h[j],h[k])+cross_y(m[i],m[j],m[k],hp[i],hp[j],hp[k]));
            jtn[k] = coeff*(cross_z(mp[i],mp[j],mp[k],h[i],h[j],h[k])+cross_z(m[i],m[j],m[k],hp[i],hp[j],hp[k]));
        }else{
            jtn[i] = 0;
            jtn[j] = 0;
            jtn[k] = 0;
        }

        double mm = m[i]*m[i] + m[j]*m[j] + m[k]*m[k];
        double mh = m[i]*h[i] + m[j]*h[j] + m[k]*h[k];
        double mhp = m[i]*hp[i] + m[j]*hp[j] + m[k]*hp[k];
        double mph = mp[i]*h[i] + mp[j]*h[j] + mp[k]*h[k];
        double mmp = m[i]*mp[i] + m[j]*mp[j] + m[k]*mp[k];

        jtn[i] += alpha[i]*coeff*((mph+mhp)*m[i]+mh*mp[i]-2*mmp*h[i]-mm*hp[i]);
        jtn[j] += alpha[i]*coeff*((mph+mhp)*m[j]+mh*mp[j]-2*mmp*h[j]-mm*hp[j]);
        jtn[k] += alpha[i]*coeff*((mph+mhp)*m[k]+mh*mp[k]-2*mmp*h[k]-mm*hp[k]);

        
        if (default_c>0){
            jtn[i] += default_c *((1-mm)*mp[i]-2*mmp*m[i]);
            jtn[j] += default_c *((1-mm)*mp[j]-2*mmp*m[j]);
            jtn[k] += default_c *((1-mm)*mp[k]-2*mmp*m[k]);
        }

    }

}

