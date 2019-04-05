#include<iostream>
extern "C" {
#include "operators.h"
}
#include "tree.hpp"
#include "utils.hpp"
#include "calculate.hpp"

void P2P(double x, double y, double z, double mux, double muy, double muz, double *F) {
  double R = sqrt(x*x + y*y + z*z);
  double R3 = R*R*R;
  double R5 = R3*R*R;
  double mu_dot_r = mux*x + muy*y + muz*z;
  F[0] += mu_dot_r / R3;
  F[1] += (3*mu_dot_r * x / R5 - mux / R3);
  F[2] += (3*mu_dot_r * y / R5 - muy / R3);
  F[3] += (3*mu_dot_r * z / R5 - muz / R3);
}

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
		  unsigned int cell, unsigned int ncrit, unsigned int exporder) {
  if (cells[cell].nleaf > ncrit) {
    for (unsigned int octant = 0; octant < 8; octant++) {
      if (cells[cell].nchild & (1 << octant)) {
	evaluate_P2M(particles, cells, cells[cell].child[octant], ncrit, exporder);
      }
    }
  }
  else {
    double *M = new double[Nterms(exporder+1)]();
    for(unsigned int i = 0; i < (cells[cell].nleaf); i++) {
      int l = cells[cell].leaf[i];
      M[1] = particles[l].mux;
      M[2] = particles[l].muy;
      M[3] = particles[l].muz;
      double dx = (particles[l].x - cells[cell].x);
      double dy = (particles[l].y - cells[cell].y);
      double dz = (particles[l].z - cells[cell].z);
      M2M(-dx, -dy, -dz,
	  M, cells[cell].M.data(),
	  exporder);
    }
    delete[] M;
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
        int p = cells[i].parent;
        double dx = cells[p].x - cells[i].x;
        double dy = cells[p].y - cells[i].y;
        double dz = cells[p].z - cells[i].z;
        M2M(dx, dy, dz, cells[i].M.data(), cells[p].M.data(), exporder);
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
                dx = particles[i].x - cells[c].x;
                dy = particles[i].y - cells[c].y;
                dz = particles[i].z - cells[c].z;
                r = sqrt(dx*dx + dy*dy + dz*dz);
                // Apply the Barnes-Hut criterion:
                if (cells[c].r > theta * r) {
                    // If the cell is 'near':
                    evaluate_M2P_and_P2P(particles, c, i, cells, F, n_crit, theta, exporder);
                }
                else {
                    // If the cell is 'far', calculate the potential and field
                    // on the particle from the multipole expansion:
		  double Fval[3] = {0.0};
		  
		  // M2P(dx, dy, dz, cells[c].M.data(), &F[4*i], exporder);
		  M2P(dx, dy, dz, cells[c].M.data(), Fval, exporder);

		  F[3*i+0] -= Fval[0];
		  F[3*i+1] -= Fval[1];
		  F[3*i+2] -= Fval[2];
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
              dx = particles[i].x - particles[j].x;
              dy = particles[i].y - particles[j].y;
              dz = particles[i].z - particles[j].z;
              // Compute the field:
              P2P(dx, dy, dz, particles[j].mux, particles[j].muy, particles[j].muz, &F[3*i]);
            }
        }
    }
}

void evaluate_approx(std::vector<Particle> &particles, std::vector<Cell> &cells, double *F, unsigned int n_crit, double theta, unsigned int exp_order) {
  #pragma omp parallel for
  for(unsigned int i = 0; i < particles.size(); i++) {
    evaluate_M2P_and_P2P(particles, 0, i, cells, F, n_crit, theta, exp_order);
  }
}


void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &F) {
    #pragma omp parallel for
    for (unsigned int i = 0; i < particles.size(); i++) {
        for (unsigned int j = 0; j < particles.size(); j++) {
            if (i != j) {
              double dx = particles[i].x - particles[j].x;
              double dy = particles[i].y - particles[j].y;
              double dz = particles[i].z - particles[j].z;
              // calculation of R and R3 will be inlined by compiler
              // so no need to worry about that.
              P2P(dx, dy, dz, particles[j].mux, particles[j].muy, particles[j].muz, &F[3*i]);
            }
        }
    }
}
