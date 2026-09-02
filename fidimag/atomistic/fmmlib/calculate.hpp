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

template <int D>
void M_sanity_check(const std::vector<Cell<D>> &cells);

template <int D>
void evaluate_P2M(std::vector<Cell<D>> &cells,
                  const double *bx, const double *by, const double *bz,
                  const double *bS, size_t ncrit, size_t exporder);

template <int D>
void evaluate_M2M(std::vector<Cell<D>> &cells,
                  const std::vector<std::vector<size_t>> &levels, size_t exporder);


template <int D>
void evaluate_L2L(std::vector<Cell<D>> &cells, const std::vector<std::vector<size_t>> &levels,
                  size_t exporder);

template <int D>
void evaluate_L2P(std::vector<Cell<D>> &cells,
                  const double *bx, const double *by, const double *bz,
                  double *F, size_t ncrit, size_t exporder);

template <int D>
void evaluate_direct(const double *bx, const double *by, const double *bz,
                     const double *bS, double *F, size_t n);


template <int D>
void interact_dehnen_lazy(const size_t A, const size_t B, const std::vector<Cell<D>> &cells,
                          const std::vector<Particle<D>> &particles,
                          const double theta, const size_t order, const size_t ncrit,
                          std::vector<Interaction> &M2L_list,
                          std::vector<Interaction> &P2P_list);




template <int D>
void evaluate_M2L_lazy(std::vector<Cell<D>> &cells,
                    std::vector<Interaction> &M2L_list,
                    std::vector<size_t> &M2L_group, size_t order);

template <int D>
void evaluate_P2P_lazy(std::vector<Cell<D>> &cells,
                    const double *bx, const double *by, const double *bz,
                    const double *body_S,
                    std::vector<Interaction> &P2P_list,
                    std::vector<size_t> &P2P_group, double *F);

template <int D>
void evaluate_M2P_and_P2P(const double *bx, const double *by, const double *bz,
  const double *bS, unsigned int p, size_t m,
  std::vector<Cell<D>> &cells, double *F, unsigned int n_crit, double theta,
  unsigned int exporder);
