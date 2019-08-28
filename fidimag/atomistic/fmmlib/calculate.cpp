#include "calculate.hpp"
#include "operators.h"
#include "tree.hpp"
#include "utils.hpp"
#include <iostream>
#include <stack>
#include <cmath>

void M_sanity_check(const std::vector<Cell> &cells) {
	double M0 = 0;
	for(size_t c = 1; c < cells.size(); c++) {
      if (cells[c].nchild == 0) {
		    M0 += cells[c].M[0];
	    }
  }
	std::cout << "Cell 0 has M0 = " << cells[0].M[0] << std::endl;
	std::cout << "Other cells   = " << M0 << std::endl;
  if (std::abs((cells[0].M[0] - M0)/M0) > 10e-10) {
    throw std::runtime_error("M0 sanity check failed");
  }
}


void P2P_Cells(size_t A, size_t B, std::vector<Cell> &cells, std::vector<Particle> &particles, double *F) {
  // A - target
  // B - source
  for (size_t p1 = 0; p1 < cells[A].nleaf; p1++) {
    //double F_p[FMMGEN_OUTPUTSIZE] = {0.0};
    size_t l1 = cells[A].leaf[p1];
    for (size_t p2 = 0; p2 < cells[B].nleaf; p2++) {
      size_t l2 = cells[B].leaf[p2];
      if (l2 != l1) {
      	double dx = (particles[l1].r[0] - particles[l2].r[0]);
      	double dy = (particles[l1].r[1] - particles[l2].r[1]);
      	double dz = (particles[l1].r[2] - particles[l2].r[2]);
      	P2P(dx, dy, dz, particles[l2].S, &F[FMMGEN_OUTPUTSIZE*l1]);
      }
    }
  }
}

void interact_dehnen(size_t A, size_t B, std::vector<Cell> &cells, std::vector<Particle> &particles,
		     double theta, size_t order, size_t ncrit, double *F) {
  double dx = (cells[A].x - cells[B].x);
  double dy = (cells[A].y - cells[B].y);
  double dz = (cells[A].z - cells[B].z);
  double R = sqrt(dx*dx + dy*dy + dz*dz);

  if (R*theta > (cells[A].rmax + cells[B].rmax)) {
    M2L(dx, dy, dz, cells[B].M, cells[A].L, order);
  }

  else if (cells[A].nchild == 0 && cells[B].nchild == 0) {
    if (cells[B].nleaf >= ncrit) {
      M2L(dx, dy, dz, cells[B].M, cells[A].L, order);
    }
    else {
      P2P_Cells(A, B, cells,particles, F);
    }
  }

  else if (cells[B].nchild == 0 || (cells[A].rmax >= cells[B].rmax && cells[A].nchild != 0)) {
      for(int oa = 0; oa < 8; oa++) {
        if (cells[A].nchild & (1 << oa)) {
          int a = cells[A].child[oa];
          interact_dehnen(a, B, cells, particles, theta, order, ncrit, F);
        }
      }
  }

  else {
    for(int ob = 0; ob < 8; ob++) {
      if (cells[B].nchild & (1 << ob)) {
        const int b = cells[B].child[ob];
        interact_dehnen(A, b, cells, particles, theta, order, ncrit, F);
      }
    }
  }
}

void interact_dehnen_lazy(const size_t A, const size_t B,
                          const std::vector<Cell> &cells,
                          const std::vector<Particle> &particles,
			                    const double theta, const size_t order,
                          const size_t ncrit,
                          std::vector<std::pair<size_t, size_t>> &M2L_list,
                          std::vector<std::pair<size_t, size_t>> &P2P_list) {
  const double dx = cells[A].x - cells[B].x;
  const double dy = cells[A].y - cells[B].y;
  const double dz = cells[A].z - cells[B].z;
  const double R = sqrt(dx*dx + dy*dy + dz*dz);

  if (R*theta > (cells[A].rmax + cells[B].rmax)) {
    //if (cells[A].nleaf < ncrit && cells[B].nleaf < ncrit) {
      std::pair<size_t, size_t> m2l_pair = std::make_pair(B, A);
      M2L_list.push_back(m2l_pair);
    //}
  }

  else if (cells[A].nchild == 0 && cells[B].nchild == 0) {
    if (cells[B].nleaf >= ncrit) {
      std::pair<size_t, size_t> m2l_pair = std::make_pair(B, A);
      M2L_list.push_back(m2l_pair);
      M2L(dx, dy, dz, cells[B].M, cells[A].L, order);
    }
    else {
      //if (cells[A].nleaf < ncrit and cells[B].nleaf < ncrit) {
    	std::pair<size_t, size_t> p2p_pair = std::make_pair(A, B);
    	P2P_list.push_back(p2p_pair);
      //}
    }
  }

  else if (cells[B].nchild == 0 || (cells[A].rmax >= cells[B].rmax && cells[A].nchild != 0)) {
    for(int oa = 0; oa < 8; oa++) {
      // For all 8 children of A, if child exists
      if (cells[A].nchild & (1 << oa)) {
    	int a = cells[A].child[oa];
    	interact_dehnen_lazy(a, B, cells, particles, theta, order, ncrit, M2L_list, P2P_list);
      }
    }
  }

  else {
    for(int ob = 0; ob < 8; ob++) {
      // for all 8 children of B, if child exists:
      if (cells[B].nchild & (1 << ob)) {
        int b = cells[B].child[ob];
        interact_dehnen_lazy(A, b, cells, particles, theta, order, ncrit, M2L_list, P2P_list);
      }
    }
  }
}

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
              size_t cell, size_t ncrit, size_t exporder) {
  #pragma omp for
  for(size_t c = 0; c < cells.size(); c++) {
    //std::cout << "Cell " << c << std::endl;
    //std::cout << "  Msize = " << Msize(exporder, FMMGEN_SOURCEORDER) << std::endl;
    size_t msize = Msize(exporder, FMMGEN_SOURCEORDER);
    double *M = new double[msize]();

    if (cells[c].nleaf < ncrit) {
      for(size_t i = 0; i < cells[c].nleaf; i++) {
        size_t l = cells[c].leaf[i];
        // Walter dehnen's definition:
        // (-1)^m / m! (x_a - z_a)^m
        double dx = (cells[c].x - particles[l].r[0]);
        double dy = (cells[c].y - particles[l].r[1]);
        double dz = (cells[c].z - particles[l].r[2]);
        for(int k = 0; k < FMMGEN_SOURCESIZE; k++) {
          // std::cout << particles[l].S[k] << std::endl;
          M[k] = particles[l].S[k];
        }
        M2M(dx, dy, dz, M, cells[c].M, exporder);
      }
   }
   delete[] M;
  }
}

void evaluate_M2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  size_t exporder) {
  /*
  evaluate_M2M(particles, cells)

  This function evaluates the multipole to
  multipole kernel. It does this by working up the
  tree from the leaf nodes, which is possible
  by iterating backwards through the nodes because
  of the way the tree is constructed.
  */
  // #pragma omp for schedule(dynamic)
  // Can't currently go up the tree in parallel.
  // Needs to be recursive or summing not correct.

  // Dehnen definition:
  // M_m(z_p) = (z_p - z_c)^n / n! M_{m - n}
  for (size_t c = cells.size() - 1; c > 0; c--) {
    size_t p = cells[c].parent;
    double dx = (cells[p].x - cells[c].x);
    double dy = (cells[p].y - cells[c].y);
    double dz = (cells[p].z - cells[c].z);
    M2M(dx, dy, dz, cells[c].M, cells[p].M, exporder);
  }
}


void evaluate_M2L_lazy(std::vector<Cell> &cells,
                       std::vector<std::pair<size_t, size_t>> &M2L_list, size_t order) {
    #pragma omp for
    for(size_t i = 0; i < M2L_list.size(); i++) {
    	size_t B = M2L_list[i].first;
    	size_t A = M2L_list[i].second;
      // Dehnen definition:
      // F_n(z_B) = M_m(z_A) * D_{n+m} (z_B - z_A)

      // So here, we've reversed the order:
      // F_n(z_A) = M_m(z_B) * D_{n+m} (z_A - z_B)
    	double dx = cells[A].x - cells[B].x;
    	double dy = cells[A].y - cells[B].y;
    	double dz = cells[A].z - cells[B].z;
    	M2L(dx, dy, dz, cells[B].M, cells[A].L, order);
    }
}

void evaluate_P2P_lazy(std::vector<Cell> &cells, std::vector<Particle> &particles,
                       std::vector<std::pair<size_t, size_t>> &P2P_list, double *F) {
   #pragma omp for
   for(size_t i = 0; i < P2P_list.size(); i++) {
       size_t A = P2P_list[i].first;
       size_t B = P2P_list[i].second;
       P2P_Cells(A, B, cells, particles, F);
   }
}


void evaluate_L2L(std::vector<Cell> &cells, size_t exporder) {
  // Can't currently go down the tree in parallel!
  // needs to be recursive or summing not correct.
  for (size_t p = 0; p < cells.size(); p++) {
    for (int octant = 0; octant < 8; octant++) {
      if (cells[p].nchild & (1 << octant)) {
        // for child c in cell p
        size_t c = cells[p].child[octant];
        double dx = cells[c].x - cells[p].x;
        double dy = cells[c].y - cells[p].y;
        double dz = cells[c].z - cells[p].z;
        L2L(dx, dy, dz, cells[p].L, cells[c].L, exporder);
      }
    }
  }
}

void evaluate_L2P(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  double *F, size_t ncrit, size_t exporder) {
  #pragma omp for schedule(runtime)
  for (size_t i = 0; i < cells.size(); i++) {
    if (cells[i].nleaf < ncrit) {
      for (size_t p = 0; p < cells[i].nleaf; p++) {
	    size_t k = cells[i].leaf[p];
        double dx = particles[k].r[0] - cells[i].x;
        double dy = particles[k].r[1] - cells[i].y;
        double dz = particles[k].r[2] - cells[i].z;
        L2P(dx, dy, dz, cells[i].L, &F[FMMGEN_OUTPUTSIZE*k], exporder);
      }
    }
  }
}

void evaluate_direct(std::vector<Particle> &particles, double *F, size_t n) {
  #pragma omp parallel for schedule(runtime)
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      if (i != j) {
      	double dx = particles[i].r[0] - particles[j].r[0];
      	double dy = particles[i].r[1] - particles[j].r[1];
      	double dz = particles[i].r[2] - particles[j].r[2];
      	P2P(dx, dy, dz, particles[j].S, &F[FMMGEN_OUTPUTSIZE*i]);
      }
    }
  }
}

void evaluate_M2P_and_P2P(std::vector<Particle> &particles, unsigned int p, unsigned int i,
  std::vector<Cell> &cells, double *F, unsigned int n_crit, double theta,
  unsigned int exporder) {
  // For particle i, start at cell p
  double dx, dy, dz, r;
  int c;
  unsigned int j;
  // If cell p is not a leaf cell:
  if (cells[p].nleaf >= n_crit) {
    // Iterate through it's children
    for (unsigned int octant = 0; octant < 8; octant++) {
      // If a child exists in a given octant:
      if (cells[p].nchild & (1 << octant)) {
        // Get the child's index
        c = cells[p].child[octant];
        // Calculate the distance from the particle to the child cell
        dx = particles[i].r[0] - cells[c].x;
        dy = particles[i].r[1] - cells[c].y;
        dz = particles[i].r[2] - cells[c].z;
        r = sqrt(dx*dx + dy*dy + dz*dz);
        // Apply the Barnes-Hut criterion:
        if (cells[c].r > theta * r) {
            // If the cell is 'near':
            evaluate_M2P_and_P2P(particles, c, i, cells, F, n_crit, theta, exporder);
        }
        else {
            // If the cell is 'far', calculate the potential and field
            // on the particle from the multipole expansion:
            M2P(dx, dy, dz, cells[c].M, &F[FMMGEN_OUTPUTSIZE*i], exporder);
        }
      }
    }
  }
  else {
    // loop in leaf cell's particles
    for(unsigned int l = 0; l < (cells[p].nleaf); l++) {
      // Get the particle index:
      j = cells[p].leaf[l];
      if (i != j) {
        // Calculate the interparticle distance:
        dx = particles[i].r[0] - particles[j].r[0];
        dy = particles[i].r[1] - particles[j].r[1];
        dz = particles[i].r[2] - particles[j].r[2];
        // Compute the field:
        P2P(dx, dy, dz, particles[j].S, &F[FMMGEN_OUTPUTSIZE*i]);
      }
    }
  }
}
