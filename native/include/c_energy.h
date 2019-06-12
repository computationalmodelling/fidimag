#pragma once
#include<cmath>
#include<iostream>
#include "c_micro_sim.h"

class Energy {
public:
    Energy() {};
    virtual ~Energy() {std::cout << "Killing A\n";};
    int interaction_id;
    bool set_up;
    int nx, ny, nz, n;
    double dx, dy, dz;
    double unit_length;
    double *spin;
    double *Ms;
    double *Ms_inv;
    double *field;
    double *energy;
    double *coordinates;
    int *ngbs;
    double compute_energy();
    // Will get the parameters from a simulation class
    void _setup(MicroSim * sim);
    virtual void compute_field(double t) {};
};


class ExchangeEnergy : public Energy {
public:
    ExchangeEnergy() {
        std::cout << "Instatiating ExchangeEnergy class; at " << this << "\n";
        this->interaction_id = 1;
    };
    ~ExchangeEnergy() {std::cout << "Killing Exchange\n";};
    void setup(double *A, MicroSim * sim) {
        _setup(sim);
        this->A = A;
    }
    double *A;
    void compute_field(double t);
};


// class AnisotropyEnergy : public Energy {
// public:
//     AnisotropyEnergy() {
//         std::cout << "Instatiating AnisotropyEnergy class; at " << this << "\n";
//         this->interaction_id = 1;
//     };
//     ~AnisotropyEnergy() {std::cout << "Killing Anisotropy\n";};
//     void init(double *Ku) {
//         this->set_up = false;
//         this->Ku = Ku;
//     }
//     double *A;
//     void compute_field(double t);
// };
