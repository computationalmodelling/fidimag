#pragma once
#include<iostream>
#include<functional>
#include "c_micro_sim.h"

enum IntegratorID {
    RK4
};

// In case we want to use conditional statements for LLG eq:
// enum LLGTerms {
//     PRECESSION,
//     DAMPING,
//     STT_ZHANGLI,
//     STT_SLONCZEWSKI,
//     STT_FIELDLIKE
// };

class Integrator;  // forward declaration
// Base class for the drivers -> just declare LLG driver for now
class MicroLLGDriver {
public:
    MicroLLGDriver() {};
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
    //
    
    // std::function<void (double * m, unsigned int n, double t)> generate_RHS_function(void);

    // TODO: declare all remaining values from the full form LLG eq
    double * alpha;
    double gamma;
    double t;
    double dt;
    MicroSim * sim;
    Integrator * integrator;

    // Will get the parameters from a simulation class
    void _setup(MicroSim * sim,
                double * alpha, double gamma, 
                // double * zl_jx, double * zl_jy, double * zl_jz, double zl_p, double zl_beta, double zl_u0,
                // double * cpp_p, double * cpp_aJ, double * cpp_beta,
                double t, double dt);
    void add_integrator(IntegratorID integrator_id);
    void run_until(double t);
    void compute_RHS(double * m, unsigned int n, double t);
};

// class MicroLLGDriver: public Driver {
// public:
//     MicroLLGDriver() {};
//     // ~MicroLLGDriver() {std::cout << "Killing base Driver\n";};
//
//
// }

// ----------------------------------------------------------------------------

// Abstract class
class Integrator {
public:
    Integrator() {};

    // We probably should better use a constructor:
    // Accepts N variables (e.g. 3 * n_spins) and a pointer to a function
    virtual void _setup(unsigned int N,
                        double dt,
                        void (*f) (double * m, unsigned int n, double t)
                        // double * init_m
                        ) {};
    // Allow magnetisation to update (via LLG + effective field)
    // Here we store a pointer to the RHS calculation function set up by the
    // driver
    void (*compute_RHS) (double * m, unsigned int n, double t) {};
    virtual void integration_step(double * m) {};

    IntegratorID integrator_id;
};

// Integrators should be independent of any driver/sim class
class Integrator_RK4: public Integrator {
public:
    Integrator_RK4() {};
    virtual ~Integrator_RK4() {std::cout << "Killing RK4 integrator\n";};

    // Will get the parameters from a simulation class
    // void _setup(MicroSim * sim, Driver * driver);
    std::vector<double> rksteps;  // N len vector, we could also use an array

    void _setup(unsigned int N,
                double dt,
                void (*f) (double * m, unsigned int n, double t)
                // double * init_m
                );
    void integration_step(double * m);

    unsigned int step_n;
    unsigned int N;
    double t;
    double dt;
};
