#include "m_driver.h"

void LLG_RHS(double * m, unsigned int n, double t, double * alpha, double gamma) {
    // Compute LLG RHS in m
}

void MicroLLGDriver::_setup(MicroSim * sim,
                            double * alpha, double gamma, double t, double dt) {

    this->sim = sim;
    this->alpha = alpha;
    this->gamma = gamma;
    this->t = t;
    this->dt = dt;

}

// TODO: Set up LLG equation function: compute_RHS
void MicroLLGDriver::compute_RHS(double * m, unsigned int n, double t) {
 
}

void MicroLLGDriver::add_integrator(IntegratorID integrator_id) {

    switch(integrator_id) {
        case RK4:
            this->integrator = new Integrator_RK4();
    }

    // Cannot use a capuring lambda (to use parameters) to make a function pointer
    // We can use a static function and create a capture-less lambda
    // See: https://deviorel.wordpress.com/2015/01/27/obtaining-function-pointers-from-lambdas-in-c/
    // Here it is ok as we do not need more than one RHS function
    // TODO: a functional approach might be better, check integrators in ANSI C
    auto _compute_RHS = [this](double * m, unsigned int n, double t) { 
        return LLG_RHS(m, n, t, this->alpha, this->gamma); 
    };
    static auto static_rhs_fun = _compute_RHS;
    void (*fptr)(double * m, 
                 unsigned int n,
                 double t) = [] (double * m, 
                                 unsigned int n,
                                 double t) {return static_rhs_fun(m, n, t);};
    this->integrator->_setup(this->sim->n, this->dt, 
                             fptr
                             // this->
                             );
}

// ----------------------------------------------------------------------------

void Integrator_RK4::_setup(unsigned int N,  // Change to unsigned long??
                            double dt,
                            // check these limits:
                            void (*f) (double * m, unsigned int n, double t)
                            // double * init_m  // Might not be necessary!
                            ) {
    // Initial values
    this->rksteps.resize(N * 4, 0.0);
    this->step_n = 0;
    this->t = 0.0;

    this->dt = dt;
    this->compute_RHS = f;

    // for (unsigned int i = 0; i < N; ++i ) {
    //     rksteps[i] = init_m[i];
    // }
}

void Integrator_RK4::integration_step(double * m) {

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
        this->compute_RHS(&rksteps[RKSTEP * N], N, t + t_factor[RKSTEP] * dt);
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
