//############################################
//#
//# Functions for running the Barnes-Hut
//# method for gravitational source particles.
//#
//# (C) Ryan Pepper, 2018
//# University of Southampton, UK
//#
//#
//###########################################
#pragma once
#include "tree.hpp"
#include "utils.hpp"
#include <iostream>

void P2P(double x, double y, double z, double mux, double muy, double muz, double *F);

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
		  size_t cell, size_t ncrit, size_t exporder);

void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  size_t exporder);


void evaluate_L2L(std::vector<Cell> &cells, size_t exporder);

void evaluate_L2P(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  double *F, size_t ncrit, size_t exporder);

void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &F);
// void evaluate_direct(std::vector<Particle> &particles, std::vector<double>
// &Bx, std::vector<double> &By, std::vector<double> &Bz);

void interact_dehnen(size_t A, size_t B, std::vector<Cell> &cells, std::vector<Particle> &particles, double theta, size_t order, size_t ncrit, double *F);
