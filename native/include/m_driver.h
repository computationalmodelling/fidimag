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
    MicroLLGDriver(MicroSim * sim, double * alpha, double gamma,
                   double t, double dt) {
        sim = sim;
        alpha = alpha;
        gamma = gamma;
        t = t;
        dt = dt;
    }

    // How to pass MEMBER FUNCTION pointer to the integrator ???
    // void compute_RHS(double * m, unsigned int n, double t);
    
    // std::function<void (double * m, unsigned int n, double t)> generate_RHS_function(void);

    // TODO: declare all remaining values from the full form LLG eq
    // No problem declaring these variables as static since we allow only a
    // single instance of the MicroLLGDriver per simulation
    // Alternatively, we could try making the compute_RHS dependant on a struct
    // containing all the LLG parameters
    // double * alpha;
    // double gamma;
    double t;
    double dt;
    MicroSim * sim;
    // Use integrator with the LLG_params struct
    Integrator<LLG_params> * integrator;
    LLG_params * llg_params;
    double * dmdt;  // N len vector, we could also use an array

    // Will get the parameters from a simulation class
    void _setup(MicroSim * sim,
                double * alpha, double gamma, 
                // double * zl_jx, double * zl_jy, double * zl_jz, double zl_p, double zl_beta, double zl_u0,
                // double * cpp_p, double * cpp_aJ, double * cpp_beta,
                double t, double dt);
    void add_integrator(IntegratorID integrator_id);
    void run_until(double t);
    void compute_RHS(double * m, double * dmdt, unsigned int n, double t);
};

// ----------------------------------------------------------------------------

// Abstract class for the integrator
// The template refers to the type of the pointer to pass extra parameters to
// the integrator, e.g. the LLG Driver uses a pointer to the LLG_params struct
template <class T>
class Integrator {
public:
    Integrator() {};

    // We should probably better use a constructor:
    // Accepts N variables (e.g. 3 * n_spins) and a pointer to a function
    virtual void _setup(unsigned int N,
                        double dt
                        ) {};
    // Allow magnetisation to update (via LLG + effective field)
    // Should we store a pointer to the update function??
    // void (*compute_RHS) (double * m, double * dmdt, unsigned int n, double t) {};

    // Templated pointer for parameters, e.g. for LLG we need eff field, alpha, etc
    virtual void integration_step(void (*f) (double * m, double * dmdt, unsigned int n, double t, void * params),
                                  double * m, double * dmdt, T * params) {};

    IntegratorID integrator_id;
};

// Integrators should be independent of any driver/sim class
template <class T>
class Integrator_RK4: public Integrator<T> {
public:
    Integrator_RK4() {};
    virtual ~Integrator_RK4() {std::cout << "Killing RK4 integrator\n";};

    std::vector<double> rksteps;  // N len vector, we could also use an array

    void _setup(unsigned int N, double dt);
    void integration_step(void (*f) (double * m, double * dmdt, unsigned int n, double t, T * params),
                          double * m, double * dmdt, T * params);

    unsigned int step_n;
    unsigned int N;
    double t;
    double dt;
};
