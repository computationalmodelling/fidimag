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
#include<omp.h>

void M_sanity_check(const std::vector<Cell> &cells);

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
		  size_t cell, size_t ncrit, size_t exporder);

void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  size_t exporder);


void evaluate_L2L(std::vector<Cell> &cells, size_t exporder);

void evaluate_L2P(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  double *F, size_t ncrit, size_t exporder);

void evaluate_direct(std::vector<Particle> &particles, double *F, size_t Nparticles);

void interact_dehnen(size_t A, size_t B, std::vector<Cell> &cells, std::vector<Particle> &particles, double theta, size_t order, size_t ncrit, double *F);

void interact_dehnen_lazy(const size_t A, const size_t B, const std::vector<Cell> &cells, const std::vector<Particle> &particles,
			  const double theta, const size_t order, const size_t ncrit,
			  std::vector<std::pair<size_t, size_t>> &M2L_list,
			  std::vector<std::pair<size_t, size_t>> &P2P_list);

void P2P_Cells(size_t A, size_t B, std::vector<Cell> &cells,
	 			       std::vector<Particle> &particles, double *F);

void evaluate_P2P_lazy(std::vector<Cell> &cells,
                      std::vector<std::pair<size_t, size_t>> &P2P_list);

void evaluate_M2L_lazy(std::vector<Cell> &cells,
                     std::vector<std::pair<size_t, size_t>> &M2L_list,
                     std::vector<omp_lock_t> &M2L_locks);

void evaluate_M2L_lazy(std::vector<Cell> &cells,
                    std::vector<std::pair<size_t, size_t>> &M2L_list, size_t order);

void evaluate_P2P_lazy(std::vector<Cell> &cells, std::vector<Particle> &particles,
                    std::vector<std::pair<size_t, size_t>> &P2P_list, double *F);

void evaluate_M2P_and_P2P(std::vector<Particle> &particles, unsigned int p, unsigned int i,
   std::vector<Cell> &cells, double *F, unsigned int n_crit, double theta,
   unsigned int exporder);
