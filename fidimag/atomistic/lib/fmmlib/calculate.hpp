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

#include<iostream>
#include "expansions.hpp"
#include "tree.hpp"
#include "utils.hpp"
#include<omp.h>

#ifndef FMMLIB_CALCULATE_HPP
#define FMMLIB_CALCULATE_HPP


void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells, unsigned int cell, unsigned int ncrit, unsigned int exporder);

void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells, unsigned int exporder);

void evaluate_M2P_and_P2P(std::vector<Particle> &particles, unsigned int p,
                          unsigned int i, std::vector<Cell> &cells, std::vector<double> &Bx,
                          std::vector<double> &By, std::vector<double> &Bz, unsigned int n_crit,
                          double theta, unsigned int exporder, std::vector<std::vector<double>> &mempool);

void evaluate_approx(std::vector<Particle> &particles, std::vector<Cell> &cells, std::vector<double> &Bx,
                     std::vector<double> &By, std::vector<double> &Bz, unsigned int n_crit, double theta,
                    unsigned int exp_order);

void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &Bx, std::vector<double> &By,
                     std::vector<double> &Bz);



#endif
