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


// Abstract class
class Integrator {
public:
    Integrator() {};
    virtual void _setup(int N) {};
    // Allow magnetisation to update (via LLG + effective field) 
    // TODO: figure out how to set the update function in Driver class
    virtual void compute_RHS(void (* f) (double * m, double t)) {};
}

// Integrators should be independent of any driver/sim class
class Integrator_RK4: public Integrator {
public:
    Integrator_RK4() {};
    virtual ~Integrator_RK4() {std::cout << "Killing RK4 integrator\n";};
    // Will get the parameters from a simulation class
    // void _setup(MicroSim * sim, Driver * driver);
    std::vector<double> rk1;  // N array
    std::vector<double> rk2;  // N array
    std::vector<double> rk3;  // N array
    std::vector<double> rk4;  // N array
    void integration_step(double (*f)(double t, double m));  // compute_RHS ??
    void _setup(int N);
    int N;
    double t;
    double dt;
};
