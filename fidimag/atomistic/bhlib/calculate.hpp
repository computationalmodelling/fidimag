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
#include<iostream>
#include "tree.hpp"
#include "utils.hpp"


void P2P(double x, double y, double z, double mux, double muy, double muz, double *F);

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells, unsigned int cell, unsigned int ncrit, unsigned int exporder);
void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells, unsigned int exporder);
void evaluate_M2P_and_P2P(std::vector<Particle> &particles, unsigned int p, unsigned int i, std::vector<Cell> &cells, std::vector<double> &F, unsigned int n_crit, double theta, unsigned int exporder);
void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &F);
// void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &Bx, std::vector<double> &By, std::vector<double> &Bz);

void evaluate_approx(std::vector<Particle> &particles, std::vector<Cell> &cells,
		     double *F, unsigned int n_crit,
		     double theta,
		     unsigned int exp_order);
