#include "calculate.hpp"
#include "tree.hpp"
#include "utils.hpp"

extern "C" {
#include "operators.h"
}
  
#include <iostream>
#include <stack>
#include <cmath>


void P2P(double x, double y, double z, double mux, double muy, double muz, double *F) {
  double R2 = x*x + y*y + z*z;
  double R3 = R2*sqrt(R2);
  double R5 = R3*R2;
  double mu_dot_r = mux*x + muy*y + muz*z;
  F[0] += mu_dot_r / R3;
  F[1] += (3*mu_dot_r * x / R5 - mux / R3);
  F[2] += (3*mu_dot_r * y / R5 - muy / R3);
  F[3] += (3*mu_dot_r * z / R5 - muz / R3);
}

void evaluate_P2M(std::vector<Particle> &particles, std::vector<Cell> &cells,
		  size_t cell, size_t ncrit, size_t exporder) {
  std::cout << "P2M(" << cell << ")" << std::endl;
  if (cells[cell].nleaf >= ncrit) {
    for (size_t octant = 0; octant < 8; octant++) {
      if (cells[cell].nchild & (1 << octant)) {
	evaluate_P2M(particles, cells, cells[cell].child[octant], ncrit, exporder);
      }
    }
  }
  else {
    double *M = new double[Nterms(exporder+1)]();
    for(size_t i = 0; i < (cells[cell].nleaf); i++) {
      int l = cells[cell].leaf[i];
      M[1] = particles[l].mux;
      M[2] = particles[l].muy;
      M[3] = particles[l].muz;
      double dx = (particles[l].x - cells[cell].x);
      double dy = (particles[l].y - cells[cell].y);
      double dz = (particles[l].z - cells[cell].z);
      M2M(-dx, -dy, -dz, M, cells[cell].M.data(), exporder);
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

  for (size_t i = cells.size() - 1; i > 0; i--) {
    size_t p = cells[i].parent;
    double dx = cells[p].x - cells[i].x;
    double dy = cells[p].y - cells[i].y;
    double dz = cells[p].z - cells[i].z;
    M2M(dx, dy, dz, cells[i].M.data(), cells[p].M.data(), exporder);
  }
}


void P2P_Cells(size_t A, size_t B, std::vector<Cell> &cells,
  std::vector<Particle> &particles, double *F) {
    // A - target
    // B - source

    // P2P for the pair of cells
  //std::cout << "    P2P_Cells("<<A<<","<<B<<")"<<std::endl;

  for (size_t p1 = 0; p1 < cells[A].nleaf; p1++) {
    size_t l1 = cells[A].leaf[p1];
    for (size_t p2 = 0; p2 < cells[B].nleaf; p2++) {
      size_t l2 = cells[B].leaf[p2];
      if (l2 != l1) {
      	double dx = particles[l1].x - particles[l2].x;
      	double dy = particles[l1].y - particles[l2].y;
      	double dz = particles[l1].z - particles[l2].z;
      	//std::cout << "      P2P("<<l1<<","<< l2 << ")" << std::endl;
      	P2P(dx, dy, dz, particles[l2].mux, particles[l2].muy, particles[l2].muz, &F[4 * l1]);
      }
    }
  }
}


int check_L(Cell &cell) {
  for(size_t i = 0; i < cell.L.size(); i++) {
    if (std::isnan(cell.L[i])) {
      return 1;
    }
  }
  return 0;

}


void interact_dehnen(size_t A, size_t B, std::vector<Cell> &cells, std::vector<Particle> &particles, double theta, size_t order, size_t ncrit, double *F) {
  //  std::cout << "interact_dehnen("<<A<<","<<B<<")"<<std::endl;
  double dx = cells[A].x - cells[B].x;
  double dy = cells[A].y - cells[B].y;
  double dz = cells[A].z - cells[B].z;
  double R = sqrt(dx*dx + dy*dy + dz*dz);

  if (R*theta > (cells[A].rmax + cells[B].rmax)) {
    M2L(dx, dy, dz, cells[B].M.data(), cells[A].L.data(), order);
    if (check_L(cells[A])) {
      std::cout << "A = " << A << " B = " << B << std::endl;
    }
  }

  else if (cells[A].nchild == 0 && cells[B].nchild == 0) {
    if (cells[B].nleaf >= ncrit) {
      M2L(dx, dy, dz, cells[B].M.data(), cells[A].L.data(), order);
    }
    else {
      P2P_Cells(A, B, cells,particles, F);
    }
  }

  else if (cells[B].nchild == 0 || (cells[A].rmax >= cells[B].rmax && cells[A].nchild != 0)) {
      for(int oa = 0; oa < 8; oa++) {
        // For all 8 children of A, if child exists
        if (cells[A].nchild & (1 << oa)) {
          int a = cells[A].child[oa];
          interact_dehnen(a, B, cells, particles, theta, order, ncrit, F);
        }
      }
  }

  else {
    for(int ob = 0; ob < 8; ob++) {
      // for all 8 children of B, if child exists:
      if (cells[B].nchild & (1 << ob)) {
        int b = cells[B].child[ob];
        interact_dehnen(A, b, cells, particles, theta, order, ncrit, F);
      }
    }
  }
}

void evaluate_L2L(std::vector<Cell> &cells, size_t exporder) {
  for (size_t i = 0; i < cells.size(); i++) {
    for (int octant = 0; octant < 8; octant++) {
      if (cells[i].nchild & (1 << octant)) {
        // for child in cell i
        size_t c = cells[i].child[octant];
        double dx = cells[c].x - cells[i].x;
        double dy = cells[c].y - cells[i].y;
        double dz = cells[c].z - cells[i].z;
        L2L(dx, dy, dz, cells[i].L.data(), cells[c].L.data(), exporder);
      }
    }
  }
}

void evaluate_L2P(std::vector<Particle> &particles, std::vector<Cell> &cells,
                  double *F, size_t ncrit, size_t exporder) {

  for (size_t i = 0; i < cells.size(); i++) {
    if (cells[i].nleaf < ncrit) {
      // std::cout << "cell " << i << " is a leaf " << std::endl;
      for (size_t p = 0; p < cells[i].nleaf; p++) {
	      size_t k = cells[i].leaf[p];
	    // std::cout << "L2P from " << i << " to particle " << k << std::endl;
        double dx = particles[k].x - cells[i].x;
        double dy = particles[k].y - cells[i].y;
        double dz = particles[k].z - cells[i].z;
	double Fv[4] = {0.0};
        L2P(dx, dy, dz, cells[i].L.data(), Fv, exporder);
    	F[4*k+0] -= Fv[0];
    	F[4*k+1] -= Fv[1];
    	F[4*k+2] -= Fv[2];
    	F[4*k+3] -= Fv[3];
      }
    }
  }
}

// void evaluate_approx(std::vector<Particle> &particles, std::vector<Cell>
// &cells,
//                      std::vector<double> &F, size_t n_crit, double
//                      theta, size_t exp_order) {
//   for (size_t i = 0; i < particles.size(); i++) {
//     evaluate_M2P_and_P2P(particles, 0, i, cells, F, n_crit, theta,
//     exp_order);
//   }
// }

void evaluate_direct(std::vector<Particle> &particles, std::vector<double> &F) {

    for (size_t i = 0; i < particles.size(); i++) {
        for (size_t j = 0; j < particles.size(); j++) {
            if (i != j) {
              double dx = particles[i].x - particles[j].x;
              double dy = particles[i].y - particles[j].y;
              double dz = particles[i].z - particles[j].z;
              // calculation of R and R3 will be inlined by compiler
              // so no need to worry about that.
              P2P(dx, dy, dz, particles[j].mux, particles[j].muy, particles[j].muz, &F[4*i]);
            }
        }
    }
}
