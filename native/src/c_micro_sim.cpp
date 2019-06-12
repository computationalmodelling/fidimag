#include "c_energy.h"
#include "c_micro_sim.h"
#include "c_constants.h"


void MicroSim::setup(int nx, int ny, int nz, double dx, double dy, double dz,
                     double unit_length, double *coordinates, int *ngbs,
                     double *spin, double *Ms, double *Ms_inv,
                     double *energy, double *field, int *pins
                     ) {

    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    this->n = nx * ny * nz;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->unit_length = unit_length;
    this->coordinates = coordinates;
    this->ngbs = ngbs;

    this->spin = spin;
    this->Ms = Ms;
    this->Ms_inv = Ms_inv;
    this->energy = energy;
    this->field = field;
    this->pins = pins;

    set_up = true;
}

void MicroSim::add_interaction(Energy * interaction_ptr) {

    interactions.push_back(interaction_ptr);
    // interactions_id.push_back(int_id);

}

void MicroSim::compute_effective_field(double t) {

    // Using indices
    std::vector<Energy *>::iterator it;
    for(it = interactions.begin(); it != interactions.end(); it++) {
        (*it)->compute_field(t);
    }
}
