#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include<memory>
// #define MAX_ENERGY_TERMS 20

class Energy;  // forward declaration; Energy is defined in c_energy.h
class MicroSim {
public:
    MicroSim() {std::cout << "Instatiating MicroSim class; at " << this << "\n";};
    virtual ~MicroSim() {std::cout << "Killing MicroSim\n";};
    bool set_up;

    // From the mesh
    int nx, ny, nz, n;
    double dx, dy, dz;
    double unit_length;
    double *coordinates;
    int *ngbs;

    // Arrays of material properties
    double *spin;
    double *Ms;
    double *Ms_inv;
    double *energy;
    double *field;
    int *pins;

    // Array with interactions
    // void * interactions;
    // int interactions_id[MAX_ENERGY_TERMS];
    // Should we use a vector of weak_ptr? : std::vector<std::weak_ptr<Energy>> interactions;
    std::vector<Energy *> interactions;

    // Not necessary at least we re creating an array of void pointers:
    // std::vector<int> interactions_id;

    // Methods
    void setup(int nx, int ny, int nz, double dx, double dy, double dz,
               double unit_length, double *coordinates, int *ngbs,
               double *spin, double *Ms, double *Ms_inv,
               double *energy, double *field, int *pins
               );

    void add_interaction(Energy * interaction);

    // void print_interactions_id() {
    //     for(int i : interactions_id) std::cout << i << "\n";
    //     for(auto i : interactions) std::cout << i << "\n";
    // }

    void compute_effective_field(double t);
};

enum EnergyTermIDs {
  NONE = 0,
  EXCHANGE,
  DMI,
  ZEEMAN,
  UNIAXIAL_ANISOTROPY
};
