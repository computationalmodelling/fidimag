#include "m_driver.h"
#include "c_vectormath.h"
#include<cmath>


// A void pointer so we can point to any struct
void LLG_RHS(double * m, double *dmdt, unsigned int n, double t,
             void * parameters) {

    LLG_params * llg_params = static_cast<LLG_params*>(parameters);
    
    // Compute LLG RHS in m
    //
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

void MicroLLGDriver::_setup(MicroSim * sim,
                            double * alpha, double gamma, double t, double dt) {

    this->sim = sim;
    this->llg_params->alpha = alpha;
    this->llg_params->gamma = gamma;
    this->t = t;
    this->dt = dt;
    this->dmdt = new double[this->sim->n];
}

// TODO: Set up LLG equation function: compute_RHS
// We use the static values of the MicroLLGDriver class
// void MicroLLGDriver::compute_RHS(double * m, double * dmdt, unsigned int n, double t) {
//     // Define LLG
//     std::cout << this->gamma << std::endl;
//     std::cout << this->alpha[0] << std::endl;
// 
//     unsigned int i, j, k;
//     double coeff, mm, mh, c;
//     double hpi, hpj, hpk;
// 
//     // #pragma omp parallel for private(i,j,k,coeff,mm, mh, c, hpi,hpj,hpk)
//     for (int id = 0; id < n; id++) {
//         // Indexes for the 3 components of the spin (magnetic moment)
//         // at the i-th lattice (mesh) site  --> x, y, z
//         i = 3 * id;
//         j = i + 1;
//         k = i + 2;
// 
//         // Pinned spins do not follow the dynamical equation
//         if (this->sim->pins[id] > 0) {
//             dmdt[i] = 0;
//             dmdt[j] = 0;
//             dmdt[k] = 0;
//             continue;
//         }
// 
//         coeff = -this->gamma / (1.0 + this->alpha[id] * this->alpha[id]);
// 
//         // Dot products
//         mm = m[i] * m[i] + m[j] * m[j] + m[k] * m[k];
//         mh = m[i] * this->sim->field[i] + m[j] * this->sim->field[j] + m[k] * this->sim->field[k];
// 
//         // Usually, m is normalised, i.e., mm=1; 
//         // so hp = mm.h - mh.m = -m x (m x h)
//         // We set here the perpendicular componenet of the field
//         // but using the (m * m) product
//         hpi = mm * this->sim->field[i] - mh * m[i];
//         hpj = mm * this->sim->field[j] - mh * m[j];
//         hpk = mm * this->sim->field[k] - mh * m[k];
// 
//         // IMPORTANT: do not ignore mm !!
//         // What we've found is that if we igonre mm, i.e. using
//         //    hpi = h[i] - mh * m[i];
//         //    hpj = h[j] - mh * m[j];
//         //    hpk = h[k] - mh * m[k];
//         // the micromagnetic standard problem 4 failed to converge (?)
//         //
//         // NOTE (Fri 08 Jul 2016 13:58): In fact, the problem converges but with 2 less
//         // decimals of accuracy, compared with the OOMMF calculation
//         double mth0 = 0, mth1 = 0, mth2 = 0;
// 
//         // The first term: m x H_eff = m x H_perp
//         // if (do_precession){
//         if (std::abs(this->gamma) > 0){  // Correct with EPS value
//             mth0 = cross_x(m[i], m[j], m[k], hpi, hpj, hpk);
//             mth1 = cross_y(m[i], m[j], m[k], hpi, hpj, hpk);
//             mth2 = cross_z(m[i], m[j], m[k], hpi, hpj, hpk);
//         }
// 
//         // The RHS term of the LLG equation
//         dmdt[i] = coeff * (mth0 - hpi * alpha[id]);
//         dmdt[j] = coeff * (mth1 - hpj * alpha[id]);
//         dmdt[k] = coeff * (mth2 - hpk * alpha[id]);
// 
//         // In future, we will try the new method to integrate the LLG equation,
//         // A mixed mid-point Runge-Kutta like scheme for the integration of
//         // Landau-Lifshitz equation Journal of Applied Physics 115, 17D101
//         // (2014) if possible, we can combine it with adaptive step size, don't
//         // know how to do but it's worth a try.
// 
//         // if (default_c < 0){
//         //     c = 6 * sqrt(dm_dt[i] * dm_dt[i] +
//         //                  dm_dt[j] * dm_dt[j] + 
//         //                  dm_dt[k] * dm_dt[k]
//         //                  );
//         // } else {
//         //     c = default_c;
//         // }
//         //printf("%0.15g   %0.15g\n", c, default_c);
// 
//         // TODO: Correct the RHS term to keep m normalised
//         c = 1.;
//         dmdt[i] += c * (1 - mm) * m[i];
//         dmdt[j] += c * (1 - mm) * m[j];
//         dmdt[k] += c * (1 - mm) * m[k];
// 
//     }
// 
// }

void MicroLLGDriver::add_integrator(IntegratorID integrator_id) {

    switch(integrator_id) {
        case RK4:
            this->integrator = new Integrator_RK4();

            // We can now call the integrator as:
            // this->integrator->integration_step(LLG_RHS,
            //                                    this->sim->spin, 
            //                                    this->dmdt,
            //                                    this->llg_params);


    }

    // Cannot use a capuring lambda (to use parameters) to make a function pointer
    // We can use a static function and create a capture-less lambda
    // See: https://deviorel.wordpress.com/2015/01/27/obtaining-function-pointers-from-lambdas-in-c/
    // Here it is ok as we do not need more than one RHS function
    // TODO: a functional approach might be better, check integrators in ANSI C
    // auto _compute_RHS = [this](double * m, unsigned int n, double t) { 
    //     return LLG_RHS(m, n, t, this->alpha, this->gamma); 
    // };
    // static auto static_rhs_fun = _compute_RHS;
    // void (*fptr)(double * m, 
    //              unsigned int n,
    //              double t) = [] (double * m, 
    //                              unsigned int n,
    //                              double t) {return static_rhs_fun(m, n, t);};
    // this->integrator->_setup(this->sim->n, this->dt, 
    //                          fptr
    //                          // this->
    //                          );

    // Alternative 2: if compute_RHS is static we can pass it but we cannot
    //                call this-> within the compute_RHS function
    // this->integrator->_setup(this->sim->n, this->dt, 
    //                          this->compute_RHS
    //                          // this->
    //                          );

    // Alternative 3: Use std::bind (replace with lambdas?)
    // using namespace std::placeholders;
    // auto _compute_rhs = std::bind(&MicroLLGDriver::compute_RHS, this, _1, _2, _3, _4);
    // this->integrator->_setup(this->sim->n, this->dt, 
    //                          _compute_rhs
    //                          // this->
    //                          );
}

// ----------------------------------------------------------------------------

void Integrator_RK4::_setup(unsigned int N,  // Change to unsigned long??
                            double dt
                            // check these limits:
                            // void (*f) (double * m, double *dmdt, unsigned int n, double t)
                            // double * init_m  // Might not be necessary!
                            ) {
    // Initial values
    this->rksteps.resize(N * 4, 0.0);
    this->step_n = 0;
    this->t = 0.0;

    this->dt = dt;
    // this->compute_RHS = f;

    // for (unsigned int i = 0; i < N; ++i ) {
    //     rksteps[i] = init_m[i];
    // }
}

void Integrator_RK4::integration_step(void (*f) (double * m,
                                                 double * dmdt,
                                                 unsigned int n,
                                                 double t,
                                                 void * params),
                                      double * m, 
                                      double * dmdt,
                                      void * params) {

    double t_factor[4] = {0.0, 0.5, 0.5, 1.0};
    double m_factor[4] = {0.0, 0.5, 0.5, 1.0};

    for (unsigned int RKSTEP = 0; RKSTEP < 4; ++RKSTEP) {

        // Re-use the rk_steps vector to apply the RK step
        // For RKSTEP = 0 we initialise the first N elements of rk_steps vector with m[i]
        // Should we separate the step calculation for every i value??
        for (unsigned int i = 0; i < this->N; ++i ) {
            unsigned int RKSTEP_idx = i + (RKSTEP) * this->N;

            rksteps[RKSTEP_idx] = m[i];

            if (RKSTEP > 0) {
                rksteps[RKSTEP_idx] += m_factor[RKSTEP] * this->dt * rksteps[i + (RKSTEP - 1) * this->N];
            }

        }
        // Update the corresponding values of the RK step vector
        f(&rksteps[RKSTEP * N], dmdt, N, t + t_factor[RKSTEP] * dt, params);
        for (unsigned int i = 0; i < this->N; ++i ) rksteps[i + RKSTEP * N] = dmdt[i];
    }

    // Update time and array
    this->t += this->dt;
    for (unsigned int i = 0; i < this->N; ++i) {
        m[i] += (1 / 6) * this->dt * (    rksteps[i] +
                                      2 * rksteps[i + N] +
                                      2 * rksteps[i + 2 * N] +
                                          rksteps[i + 3 * N]);
    }

    // // RK4 step 1
    // this->compute_RHS(rk_steps.begin(), rk_steps.begin() + N, rk_steps, t);

    // // RK4 step 2
    // for (unsigned int i = 0; i < N; ++i ) {
    //     rk_steps[i + N] = m[i] + 0.5 * this->dt * rk_steps[i];
    // }
    // this->compute_RHS(rk_steps.begin() + N, rk_steps.begin() + 2 * N,
    //                   rk_steps, t + 0.5 * dt);

    // // RK4 step 3
    // for (unsigned int i = 0; i < N; ++i ) {
    //     rk_steps[i + 2 * N] = m[i] + 0.5 * this->dt * rk_steps[N + i];
    // }
    // this->compute_RHS(rk_steps.begin() + 2 * N, rk_steps.begin() + 3 * N,
    //                   rk_steps, t + 0.5 * dt);

    // // RK4 step 4
    // for (unsigned int i = 0; i < N; ++i ) {
    //     rk_steps[i + 2 * N] = m[i] + 0.5 * this->dt * rk_steps[2 * N + i];
    // }
    // this->compute_RHS(rk_steps.begin() + 3 * N, rk_steps.begin() + 4 * N,
    //                   rk_steps, t + 1.0 * dt);

}
