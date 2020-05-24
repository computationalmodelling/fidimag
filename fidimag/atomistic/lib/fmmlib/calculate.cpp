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
#include<omp.h>
#include "calculate.hpp"


void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells, unsigned int cell, unsigned int ncrit, unsigned int exporder) {
    /*
    evaluate_P2M(particles, cells, cell, ncrit)

    This function evaluates the particle to multipole
    kernel part of the multipolar Barnes-Hut method,
    by iterating down the tree recursively, with the P2M
    method called if the curent cell is a leaf cell
    (i.e. it has no child cells) and the evaluate
    function agian otherwise
    */
    if (cells[cell].nleaf > ncrit) {
        for (unsigned int octant = 0; octant < 8; octant++) {
            if (cells[cell].nchild & (1 << octant)) {
                evaluate_P2M(particles, cells, cells[cell].child[octant], ncrit, exporder);
            }
        }
    }
    else {
        P2M(particles, cell, cells, ncrit, exporder);
    }
}


void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells, unsigned int exporder) {
    /*
    evaluate_M2M(particles, cells)

    This function evaluates the multipole to
    multipole kernel. It does this by working up the
    tree from the leaf nodes, which is possible
    by iterating backwards through the nodes because
    of the way the tree is constructed.
    */
    for(int i = cells.size() - 1; i > 0; i--) {
        M2M(cells[i].parent, i, cells, exporder);
    }
}

void evaluate_M2P_and_P2P(std::vector<Particle> &particles, unsigned int p, unsigned int i, std::vector<Cell> &cells, std::vector<double> &Bx, std::vector<double> &By, std::vector<double> &Bz, unsigned int n_crit, double theta, unsigned int exporder, std::vector<std::vector<double>> &mempool) {
    double dx, dy, dz, r;
    int c, j;
    if (cells[p].nleaf >= n_crit) {
        for (unsigned int octant = 0; octant < 8; octant++) {
            if (cells[p].nchild & (1 << octant)) {
                c = cells[p].child[octant];
                dx = particles[i].x - cells[c].x;
                dy = particles[i].y - cells[c].y;
                dz = particles[i].z - cells[c].z;
                r = sqrt(dx*dx + dy*dy + dz*dz);
                if (cells[c].r > theta * r) {
                    evaluate_M2P_and_P2P(particles, c, i, cells, Bx, By, Bz, n_crit, theta, exporder, mempool);
                }
                else {
		    int thread = omp_get_thread_num();
                    M2P(c, i, cells, particles, Bx, By, Bz, exporder, mempool[thread]);
                }
            }
        }
    }
    else {
        // loop in leaf cell's particles
        for(unsigned int l = 0; l < (cells[p].nleaf); l++) {
            j = cells[p].leaf[l];
            P2P(i, j, particles, Bx, By, Bz);
        }
    }
}

void evaluate_approx(std::vector<Particle> &particles, std::vector<Cell> &cells, std::vector<double> &Bx, std::vector<double> &By, std::vector<double> &Bz, unsigned int n_crit, double theta, unsigned int exp_order) {
  int nthreads;
  #pragma omp parallel
  {	
	#pragma omp single
	nthreads = omp_get_num_threads();
  }
  std::vector<std::vector<double>> mempool;
  for(int i = 0; i < nthreads; i++) {
	mempool.push_back(std::vector<double>(Nterms(exp_order + 2)));
  }
  #pragma omp parallel for schedule(guided)
  for(unsigned int i = 0; i < particles.size(); i++) {
    evaluate_M2P_and_P2P(particles, 0, i, cells, Bx, By, Bz, n_crit, theta, exp_order, mempool);
  }
}

void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &Bx, std::vector<double> &By, std::vector<double> &Bz) {
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int i = 0; i < particles.size(); i++) {
        for (unsigned int j = 0; j < particles.size(); j++) {
            P2P(i, j, particles, Bx, By, Bz);
        }
    }
}


