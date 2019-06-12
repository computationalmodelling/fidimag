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
    double *A;
    void setup(double *A, MicroSim * sim) {
        _setup(sim);
        this->A = A;
    }
    void compute_field(double t);
};


class AnisotropyEnergy : public Energy {
public:
    AnisotropyEnergy() {
        std::cout << "Instatiating AnisotropyEnergy class; at " << this << "\n";
        this->interaction_id = UNIAXIAL_ANISOTROPY;
    };
    ~AnisotropyEnergy() {std::cout << "Killing Anisotropy\n";};
    double *Ku;
    double *axis;
    void setup(double *Ku, double *axis, MicroSim * sim) {
        _setup(sim);
        this->Ku = Ku;
        this->axis = axis;
    }
    void compute_field(double t);
};


class DMIEnergy : public Energy {
public:
    DMIEnergy() {
        std::cout << "Instatiating DMIEnergy class; at " << this << "\n";
        this->interaction_id = DMI;
    };
    ~DMIEnergy() {std::cout << "Killing DMI\n";};
    double *D;
    double *dmi_vector;
    int n_dmis;
    void setup(double * D, double *dmi_vector, int n_dmis, MicroSim * sim) {
        _setup(sim);
        this->D = D;
        this->dmi_vector = dmi_vector;
        this->n_dmis = n_dmis;
    }
    void compute_field(double t);
};
