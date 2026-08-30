#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include<array>
#include "utils.hpp"
#include <cstdint>

/*! \brief A single cell-cell interaction.

    Replaces std::pair<size_t,size_t>, which stored (source, target) for the
    M2L list but (target, source) for the P2P list. Both lists are grouped by
    target, so having the key in a different position in each was a latent
    race: using the wrong member compiles cleanly and yields groups in which a
    target appears more than once, at which point two threads write the same
    accumulator unsynchronised. Named members make that mistake impossible.

    uint32_t is ample for cell indices -- build_tree rejects anything larger --
    and halves the memory the interaction lists occupy. */
struct Interaction {
  uint32_t target;
  uint32_t source;
};


/*! \brief Particle class used to store position and source strength. */
class Particle {
public:
  double *r; // Address of 3-vector position of the particle
  double *S; // Address of q
};

class Cell {
public:
  size_t nleaf; /*!< \brief Number of particles held in cell.
                           This counter is incremented every time a particle is
                           added to it in the
                           \ref build_tree function. This continues to be the case
                           even when the cell has been split, as we use it to keep
                           track of whether a cell has been split or not to
                           save on memory, rather than having another variable. */
  size_t nchild; /*!< \brief Number of child cells occupied.

                            Binary counter showing whether a given octant is
                            occupied by a child cell.<br>I.e. if 0001001, then there
                            are two child cells held by this cell. */
  size_t level; /*!< \brief Level of the tree that the cell sits at.

                           This is 0 for the root cell, 1 for the 1st level, etc.
                           */
  std::vector<size_t> child; /*!< \brief Indices of child octants. */
  //std::vector<double> M;
  //std::vector<double> L;
  double *M;
  double *L;
  std::vector<size_t> leaf; /*!< \brief Indices of particles in the cell. */
  /*! \brief Start of this cell's particles in the tree's Morton-ordered
      arrays. A leaf's particles occupy [body_offset, body_offset + nleaf).
      Meaningful for leaf cells only. */
  size_t body_offset;
  double x; /*!< \brief x coordinates of cell centre. */
  double y; /*!< \brief y coordinates of cell centre. */
  double z; /*!< \brief z coordinates of cell centre. */
  double r; /*!< \brief Radius of cell
                Must be sufficiently large for the root cell to bound the
                particles.

                Note: I may change this in future so it is calculated
                in build_tree rather than user specified.
                */
  double rmax;
  size_t parent; /*!< \brief Index of parent cell of this cell. */
  Cell(double x, double y, double z, double r, size_t parent, size_t level, size_t ncrit);
  ~Cell();
  Cell(const Cell& other);
  Cell(Cell&& other);
  void clear();
  void resize(size_t order);
  /*! Copy operator for the Cell class */
  Cell& operator=(const Cell& other) {
    this->nleaf = other.nleaf;
    this->nchild = other.nchild;
    this->level = other.level;
    this->child = other.child;
    this->M = other.M;
    this->L = other.L;

    this->leaf = other.leaf;
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->r = other.r;
    this->parent = other.parent;
    return *this;
  }
};


class Tree {
public:
  size_t order;
  size_t ncrit;
  double theta;
  std::vector<Particle> particles;
  std::vector<Cell> cells;
  std::vector<double> M;
  std::vector<double> L;
  std::vector<Interaction> M2L_list;
  /*! \brief Offsets into M2L_list, indexed by target cell.
      Entries whose target is cell A occupy [M2L_group[A], M2L_group[A+1]).
      Produced directly by the counting sort in build_tree, and relied on by
      evaluate_M2L_lazy to give each thread exclusive ownership of one target's
      local expansion. */
  std::vector<size_t> M2L_group;
  std::vector<Interaction> P2P_list;
  /*! \brief Offsets into P2P_list, indexed by target cell. Same role as
      M2L_group: gives each thread sole ownership of one target's output. */
  std::vector<size_t> P2P_group;

  /*! \brief Particles permuted into tree (Morton) order, stored by value.

      An octree subdivided by (x>cx) + ((y>cy)<<1) + ((z>cz)<<2) visits cells in
      Z-order, so emitting each leaf's particles in cell order makes every leaf
      a contiguous range -- no key computation needed.

      Coordinates are held BY VALUE and split into separate x/y/z arrays. The
      Particle class stores double* into caller memory, so the P2P inner loop
      previously went leaf[p] -> index -> pointer -> scattered load: three
      dependent loads per source particle. SoA gives the vectorised kernel
      unit-stride loads instead of stride-3 gathers. */
  std::vector<double> body_x, body_y, body_z;
  std::vector<double> body_S;    // FMMGEN_SOURCESIZE * nparticles
  std::vector<size_t> body_perm; // body_perm[sorted] = original index

  /*! \brief Morton-ordered output accumulator, persistent across solves.

      Was a local `std::vector<double>` inside compute_field_fmm, so every call
      paid a 32 MB allocation (at N=1M), the kernel's page faults on first
      touch, and a serial zero -- per solve, for a tree whose geometry never
      changes. The intended usage is build once and solve many times, so this is
      allocated once in build_tree and re-zeroed in parallel per solve.

      The zeroing loop uses schedule(static), which matters for more than
      overhead: the same thread re-touches the same pages on every solve, so the
      pages stay on that thread's NUMA node and its slice stays cache-warm from
      one solve to the next. */
  std::vector<double> Fm;
  std::vector<std::vector<size_t>> levels;  // Cells grouped by tree level for parallel traversal
  void compute_field_fmm(double *F);
  void compute_field_bh(double *F);
  void compute_field_exact(double *F);
  //! Scatter a Morton-ordered result into the caller's particle ordering.
  void scatter_output(const double *Fm, double *F);
  //! Re-read source values through the (live) Particle pointers into
  //! body_S. Positions and tree topology never change after build_tree, but
  //! callers that reuse one Tree across many solves -- e.g. a spin dynamics
  //! integrator calling compute_field_* once per timestep with the same
  //! geometry but a new moment array each time -- need every solve to see
  //! the current source values, not the ones snapshotted at build time.
  void refresh_sources();
private:
  void clear_M();
  void clear_L();
  void clear_Fm();
};

void printTreeParticles(std::vector<Cell> &cells, size_t cell, size_t depth);

void add_child(std::vector<Cell> &cells, int octant, size_t p, size_t ncrit);

void split_cell(std::vector<Cell> &cells, std::vector<Particle> &particles, size_t p, size_t ncrit);

Tree build_tree(double *pos, double *mu, size_t nparticles, size_t ncrit, size_t order, double theta);
