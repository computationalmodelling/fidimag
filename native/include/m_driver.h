#pragma once
#include<iostream>
#include "c_micro_sim.h"

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
};


class Integrator_RK4 {
public:
    Integrator_RK4() {};
    virtual ~Integrator_RK4() {std::cout << "Killing RK4 integrator\n";};
    // Will get the parameters from a simulation class
    // void _setup(MicroSim * sim, Driver * driver);
    std::vector<double> rk_steps;  // N * 4 array
    void integration_step(double (*f)(double t, double y));
};
