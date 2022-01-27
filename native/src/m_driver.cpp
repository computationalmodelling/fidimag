#include "m_driver.h"
#include "c_vectormath.h"
#include<cmath>


// The RHS function to integrate the LLG equation
// Extra parameters are passed to the integrator in a struct
void LLG_RHS(double * m, std::vector<double>& dmdt, unsigned int n, double t,
             LLG_params * llg_params) {

    // Only if using void pointer for parameters:
    // LLG_params * llg_params = static_cast<LLG_params*>(parameters);

    // Compute LLG RHS in m
    unsigned int i, j, k;
    double coeff, mm, mh, c;
    double hpi, hpj, hpk;

    // #pragma omp parallel for private(i,j,k,coeff,mm, mh, c, hpi,hpj,hpk)
    for (int id = 0; id < n; id++) {
        // Indexes for the 3 components of the spin (magnetic moment)
        // at the i-th lattice (mesh) site  --> x, y, z
        i = 3 * id;
        j = i + 1;
        k = i + 2;

        // Pinned spins do not follow the dynamical equation
        if (llg_params->pins[id] > 0) {
            dmdt[i] = 0;
            dmdt[j] = 0;
            dmdt[k] = 0;
            continue;
        }

        coeff = -llg_params->gamma / (1.0 + llg_params->alpha[id] * llg_params->alpha[id]);

        // Dot products
        mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
        mh = m[i] * llg_params->field[i] +
             m[j] * llg_params->field[j] +
             m[k] * llg_params->field[k];

        // Usually, m is normalised, i.e., mm=1;
        // so hp = mm.h - mh.m = -m x (m x h)
        // We set here the perpendicular componenet of the field
        // but using the (m * m) product
        hpi = mm * llg_params->field[i] - mh * m[i];
        hpj = mm * llg_params->field[j] - mh * m[j];
        hpk = mm * llg_params->field[k] - mh * m[k];

        // IMPORTAthis->NT: do not ignore mm !!
        // What we've found is that if we igonre mm, i.e. using
        //    hpi = h[i] - mh * m[i];
        //    hpj = h[j] - mh * m[j];
        //    hpk = h[k] - mh * m[k];
        // the micromagnetic standard problem 4 failed to converge (?)
        //
        // this->NOTE (Fri 08 Jul 2016 13:58): In fact, the problem converges but with 2 less
        // decimals of accuracy, compared with the OOMMF calculation
        double mth0 = 0, mth1 = 0, mth2 = 0;

        // The first term: m x H_eff = m x H_perp
        // if (do_precession){
        if (std::abs(llg_params->gamma) > 0){  // Correct with EPS value
            mth0 = cross_x(m[i], m[j], m[k], hpi, hpj, hpk);
            mth1 = cross_y(m[i], m[j], m[k], hpi, hpj, hpk);
            mth2 = cross_z(m[i], m[j], m[k], hpi, hpj, hpk);
        }

        // The RHS term of the LLG equation
        dmdt[i] = coeff * (mth0 - hpi * llg_params->alpha[id]);
        dmdt[j] = coeff * (mth1 - hpj * llg_params->alpha[id]);
        dmdt[k] = coeff * (mth2 - hpk * llg_params->alpha[id]);

        // In future, we will try the new method to integrate the LLG equation,
        // A mixed mid-point Runge-Kutta like scheme for the integration of
        // Landau-Lifshitz equation Journal of Applied Physics 115, 17D101
        // (2014) if possible, we can combine it with adaptive step size, don't
        // know how to do but it's worth a try.

        // if (default_c < 0){
        //     c = 6 * sqrt(dm_dt[i] * dm_dt[i] +
        //                  dm_dt[j] * dm_dt[j] +
        //                  dm_dt[k] * dm_dt[k]
        //                  );
        // } else {
        //     c = default_c;
        // }
        //printf("%0.15g   %0.15g\n", c, default_c);

        // TODO: Correct the RHS term to keep m normalised
        c = 1.;
        dmdt[i] += c * (1 - mm) * m[i];
        dmdt[j] += c * (1 - mm) * m[j];
        dmdt[k] += c * (1 - mm) * m[k];

    }

}

// Passes the required parameters for the integrator to the LLG_params struct
void MicroLLGDriver::setup(MicroSim * sim,
                           double * alpha, double gamma, double t, double dt) {

    this->llg_params = new LLG_params;  // initialize struct
    this->sim = sim;
    this->llg_params->gamma = gamma;
    this->llg_params->alpha = alpha;
    this->llg_params->field = sim->field;
    this->llg_params->pins = sim->pins;
    this->t = t;    // Should be handled by integrator
    this->dt = dt;  // should be handled by integrator
    this->dmdt.resize(3 * this->sim->n);

    // DEBUG
    std::cout << "[driver] Gamma:" << this->llg_params->gamma << std::endl;
    std::cout << "[driver] Sim n:" << this->sim->n << std::endl;
}


void MicroLLGDriver::add_integrator(IntegratorID integrator_id) {

    switch(integrator_id) {
        case RK4:
            // How to call Derived class integrator if the pointer is to the Base class?
            this->integrator = new IntegratorRK4<LLG_params>(this->sim->n);
            // Best create a set_integrator_parameters function
            this->integrator->setup(this->sim->n, this->dt);
            // TODO: separate function to set time?
            this->integrator->t = 0.0;
            this->integrator->step_n = 0;

            std::cout << "Setup integrator" << std::endl;

    }

}

void MicroLLGDriver::run_until(double t_final){

    while(this->integrator->t <= t_final) {
        this->integrator->integration_step(LLG_RHS,
                                           this->sim->spin,
                                           this->dmdt,
                                           this->llg_params);
    }

}

// ----------------------------------------------------------------------------

template <class T>
void IntegratorRK4<T>::setup(unsigned int N, double dt) {
    // Initial values
    this->step_n = 0;
    this->t = 0.0;
    this->dt = dt;
    this->N = N;
    // this->compute_RHS = f;

    // Initial values??
    // for (unsigned int i = 0; i < this->N; ++i ) {
    //     rksteps[i] = init_m[i];
    // }
}

template <class T>
void IntegratorRK4<T>::integration_step(void (*f) (double * m,
                                                   std::vector<double>& dmdt,
                                                   unsigned int n,
                                                   double t,
                                                   T * params),
                                        double * m,
                                        std::vector<double>& dmdt,
                                        T * params) {

    // To be defined in constructor/setup:
    double t_factor[4] = {0.0, 0.5, 0.5, 1.0};
    double m_factor[4] = {0.0, 0.5, 0.5, 1.0};
    std::vector<double>& rksteps = this->integratorData;

    for (unsigned int RKSTEP = 0; RKSTEP < 4; ++RKSTEP) {

        // Re-use the rk_steps vector to apply the RK step
        // For RKSTEP = 0 we initialise the first this->N elements of rk_steps vector with m[i]
        // Should we separate the step calculation for every i value??
        for (unsigned int i = 0; i < this->N; ++i ) {
            unsigned int RKSTEP_idx = i + (RKSTEP) * this->N;

            rksteps[RKSTEP_idx] = m[i];

            if (RKSTEP > 0) {
                rksteps[RKSTEP_idx] += m_factor[RKSTEP] * this->dt * rksteps[i + (RKSTEP - 1) * this->N];
            }

        }
        // Update the corresponding values of the RK step vector
        f(&rksteps[RKSTEP * this->N], dmdt, this->N, this->t + t_factor[RKSTEP] * this->dt, params);
        // Copy RHS values to the corresponding RK step array section
        for (unsigned int i = 0; i < this->N; ++i ) rksteps[i + RKSTEP * this->N] = dmdt[i];
    }

    // Update time and array
    this->t += this->dt;
    for (unsigned int i = 0; i < this->N; ++i) {
        m[i] += (1 / 6) * this->dt * (    rksteps[i] +
                                      2 * rksteps[i + this->N] +
                                      2 * rksteps[i + 2 * this->N] +
                                          rksteps[i + 3 * this->N]);
    }

    // // RK4 step 1
    // this->compute_RHS(rk_steps.begin(), rk_steps.begin() + this->N, rk_steps, t);

    // // RK4 step 2
    // for (unsigned int i = 0; i < this->N; ++i ) {
    //     rk_steps[i + this->N] = m[i] + 0.5 * this->dt * rk_steps[i];
    // }
    // this->compute_RHS(rk_steps.begin() + this->N, rk_steps.begin() + 2 * this->N,
    //                   rk_steps, t + 0.5 * dt);

    // // RK4 step 3
    // for (unsigned int i = 0; i < this->N; ++i ) {
    //     rk_steps[i + 2 * this->N] = m[i] + 0.5 * this->dt * rk_steps[this->N + i];
    // }
    // this->compute_RHS(rk_steps.begin() + 2 * this->N, rk_steps.begin() + 3 * this->N,
    //                   rk_steps, t + 0.5 * dt);

    // // RK4 step 4
    // for (unsigned int i = 0; i < this->N; ++i ) {
    //     rk_steps[i + 2 * this->N] = m[i] + 0.5 * this->dt * rk_steps[2 * this->N + i];
    // }
    // this->compute_RHS(rk_steps.begin() + 3 * this->N, rk_steps.begin() + 4 * this->N,
    //                   rk_steps, t + 1.0 * dt);

}

// For testing: ---------------------------------------------------------------

template class IntegratorRK4<TestParams>;

// TODO: If we move the integrator to a separate file, we should declare all
// the templates used by Fidimag here, e.g.
// template class IntegratorRK4<LLG_params>;
