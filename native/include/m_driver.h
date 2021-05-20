#pragma once
#include<iostream>
#include "c_micro_sim.h"

class Integrator;  // forward declaration
// Base class for the drivers
class Driver {
public:
    Driver() {};
    virtual ~Driver() {std::cout << "Killing base Driver\n";};
    // Will get the parameters from a simulation class
    void _setup(MicroSim * sim);

    double * alpha;
    double gamma;
    double t;
    MicroSim * sim;
    Integrator * integrator;

    virtual void run_until(double t) {};
};

class MicroLLGDriver: public Driver {
public:
    MicroLLGDriver() {};
}

// ----------------------------------------------------------------------------

// Abstract class
class Integrator {
public:
    Integrator() {};
    virtual void _setup(unsigned int N,
                        double dt,
                        void (*f) (std::vector<int>::iterator begin,
                                   std::vector<int>::iterator end,
                                   double * m, double t),
                        double * init_m) {};
    // Allow magnetisation to update (via LLG + effective field)
    // Here we store a pointer to the RHS calculation function set up by the
    // driver
    void (*compute_RHS) (std::vector<int>::iterator begin,
                         std::vector<int>::iterator end,
                         double * m, double t) {};
}

// Integrators should be independent of any driver/sim class
class Integrator_RK4: public Integrator {
public:
    Integrator_RK4() {};
    virtual ~Integrator_RK4() {std::cout << "Killing RK4 integrator\n";};
    // Will get the parameters from a simulation class
    // void _setup(MicroSim * sim, Driver * driver);
    std::vector<double> rksteps;  // N array
    void integration_step(double * m);  // compute_RHS ??
    // Accepts N variables (e.g. 3 * n_spins) and a pointer to a function
    void _setup(unsigned int N,
                double dt,
                void (*f) (std::vector<int>::iterator begin,
                           std::vector<int>::iterator end,
                           double * m, double t),
                double * init_m);
    unsigned int step_n;
    unsigned int N;
    double t;
    double dt;
};
