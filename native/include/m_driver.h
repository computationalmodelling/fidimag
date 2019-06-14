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
};
