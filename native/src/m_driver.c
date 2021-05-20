#include "m_driver.h"

void Integrator_RK4::_setup(unsigned int N,  // Change to unsigned long??
                            double dt,
                            // check these limits:
                            void (*f) (std::vector<int>::iterator begin,
                                       std::vector<int>::iterator end,
                                       double * m, double t),
                            double * init_m  // Might not be necessary!
                            ) {
    // Initial values
    this->rk_steps(N * 4, 0.0);
    this->step_n = 0;
    this->t = 0.0;

    this->dt = dt;
    this->compute_RHS = &f;

    for (unsigned int i = 0; i < N; ++i ) {
        rk_steps[i] = init_m[i];
    }
}

void Integrator_RK4::integration_step(double * m) {

    this->rk_steps(N * 4, 0.0);
    double t_factor[4] = {0.0, 0.5, 0.5, 1.0}
    double m_factor[4] = {0.0, 0.5, 0.5, 1.0}

    for (unsigned int RKSTEP = 0; RKSTEP < 4; ++RKSTEP) {

        // Re-use the rk_steps vector to apply the RK step
        // For RKSTEP = 0 we initialise the first N elements of rk_steps vector with m[i]
        // Should we separate the step calculation for every i value??
        for (unsigned int i = 0; i < this->N; ++i ) {
            unsigned int RKSTEP_idx = i + (RKSTEP) * this->N;

            rk_steps[RKSTEP_idx] = m[i];

            if (RKSTEP > 0) {
                rk_steps[RKSTEP_idx] += m_factor[RKSTEP] * this->dt * rk_steps[i + (RKSTEP - 1) * this->N];
            }

        }
        // Update the corresponding values of the RK step vector
        this->compute_RHS(rk_steps.begin() + RKSTEP * N, 
                          rk_steps.begin() + (RKSTEP + 1) * N, 
                          rk_steps, 
                          t + t_factor[RKSTEP] * dt);
    }

    this->t += this->dt;
    for (unsigned int i = 0; i < this->N; ++i) {
        m[i] += (1 / 6) * this->dt * (    rk_steps[i] + 
                                      2 * rk_steps[i + N] + 
                                      2 * rk_steps[i + 2 * N] + 
                                          rk_steps[i + 3 * N]);
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
