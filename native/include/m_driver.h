#pragma once
#include<iostream>
#include<functional>
#include "c_micro_sim.h"

enum IntegratorID {
    RK4
};

// We will set the extra parameters for the integrator as a struct
// The integrator allows m, dmdt, time, dt and extra parameters
struct LLG_params {
    double gamma;
    double * alpha;
    double * field;
    int * pins;
};

// In case we want to use conditional statements for LLG eq:
// enum LLGTerms {
//     PRECESSION,
//     DAMPING,
//     STT_ZHANGLI,
//     STT_SLONCZEWSKI,
//     STT_FIELDLIKE
// };

template <class T>
class Integrator;  // forward declaration

// Base class for the drivers -> just declare LLG driver for now
class MicroLLGDriver {
public:
    MicroLLGDriver() {std::cout << "Instatiating LLG Driver" << std::endl;};
    ~MicroLLGDriver() {std::cout << "Killing LLG Driver\n";};

    // Testing constructor
    // MicroLLGDriver(MicroSim * sim, double * alpha, double gamma,
    //                double t, double dt) {
    //     sim = sim;
    //     alpha = alpha;
    //     gamma = gamma;
    //     t = t;
    //     dt = dt;
    // }

    // TODO: declare all remaining values from the full form LLG eq
    // Parameters are set in the LLG_params struct
    LLG_params * llg_params;
    double t;
    double dt;
    MicroSim * sim;
    // This creates a pointer to the base class Integrator; with add_integrator
    // the virtual functions are override by corresp. derived class functions
    // of the chosen integrator
    Integrator<LLG_params> * integrator;
    std::vector<double> dmdt;  // N len vector, we could also use an array

    // Will get the parameters from a simulation class
    void _setup(MicroSim * sim,
                double * alpha, double gamma, 
                // double * zl_jx, double * zl_jy, double * zl_jz, double zl_p, double zl_beta, double zl_u0,
                // double * cpp_p, double * cpp_aJ, double * cpp_beta,
                double t, double dt);
    void add_integrator(IntegratorID integrator_id);
    void run_until(double t_final);
    void single_integrator_step();

    // void compute_RHS(double * m, std::vector<double>& dmdt,
    //                  unsigned int n, double t);
};

// ----------------------------------------------------------------------------

// Abstract class for the integrator
// The template refers to the type of the pointer to pass extra parameters to
// the integrator, e.g. the LLG Driver uses a pointer to the LLG_params struct
template <class T>
class Integrator {
public:
    // Integrator() = {};
    // We should probably better use a constructor:
    // Accepts N variables (e.g. 3 * n_spins) and a pointer to a function
    virtual void _setup(unsigned int N, double dt) = 0;

    // Templated pointer for parameters, e.g. for LLG we need eff field, alpha, etc
    virtual void integration_step(void (*f) (double * m, std::vector<double>& dmdt, 
                                             unsigned int n, double t, T * params),
                                  double * m, std::vector<double>& dmdt, 
                                  T * params) = 0;

    IntegratorID integrator_id;

    unsigned int step_n;
    unsigned int N;
    double t;
    double dt;

    std::vector<double> integratorData;
};

// Integrators should be independent of any driver/sim class
template <class T>
class IntegratorRK4: public Integrator<T> {
public:
    IntegratorRK4(size_t N) { 
        // Know we need 4xN for RK4
        this->integratorData.resize(4 * N);
        std::cout << "Instatiating RK4" << std::endl;
    };

    virtual ~IntegratorRK4() {std::cout << "Killing RK4 integrator\n";};

    // std::vector<double> rksteps;  // N len vector, we could also use an array

    void _setup(unsigned int N, double dt);
    void integration_step(void (*f) (double * m, std::vector<double>& dmdt, 
                                     unsigned int N, double t, T * params),
                          double * m, std::vector<double>& dmdt, T * params);
};
