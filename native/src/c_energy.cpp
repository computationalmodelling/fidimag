#include "c_energy.h"
#include "c_constants.h"

double Energy::compute_energy() {
    double sum = 0;
    for(int i = 0; i < n; i++) {
        sum += energy[i];
    }
    return sum * (dx * dy * dz * std::pow(unit_length, 3));
}

void Energy::_setup(MicroSim * sim) {

    this->nx = sim->nx;
    this->ny = sim->ny;
    this->nz = sim->nz;
    this->n = sim->nx * ny * nz;
    this->dx = sim->dx;
    this->dy = sim->dy;
    this->dz = sim->dz;
    this->unit_length = sim->unit_length;

    // Arrays
    this->spin = sim->spin;
    this->Ms = sim->Ms;
    this->Ms_inv = sim->Ms_inv;
    this->coordinates = sim->coordinates;
    this->ngbs = sim->ngbs;
    this->energy = sim->energy;
    this->field = sim->field;
    // this->pins = sim->pins;

    set_up = true;
}
