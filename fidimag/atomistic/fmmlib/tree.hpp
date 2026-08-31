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
    and halves the memory the interaction lists occupy. Dimension-independent:
    it names cells, not coordinates. */
struct Interaction {
  uint32_t target;
  uint32_t source;
};


/*! \brief Particle class used to store position and source strength.

    r points to a D-wide position (caller-owned memory); S to the
    FMMGEN_SOURCESIZE-wide source strength, which is always the full physical
    moment vector regardless of D -- a dipole in a 2D (planar) tree still has
    3 moment components, only its position is 2-wide. */
template <int D>
class Particle {
  static_assert(D == 2 || D == 3, "Tree/Cell/Particle only support D = 2 or 3");
public:
  double *r;
  double *S;
};

template <int D>
class Cell {
public:
  size_t nleaf; /*!< \brief Number of particles held in cell.
                           This counter is incremented every time a particle is
                           added to it in the
                           \ref build_tree function. This continues to be the case
                           even when the cell has been split, as we use it to keep
                           track of whether a cell has been split or not to
                           save on memory, rather than having another variable. */
  size_t nchild; /*!< \brief Bitmask of occupied child octants/quadrants.

                            D=3: up to 8 bits (octree); D=2: up to 4 bits
                            (quadtree). Bit `k` set means child `k` exists. */
  size_t level; /*!< \brief Level of the tree that the cell sits at.

                           This is 0 for the root cell, 1 for the 1st level, etc.
                           */
  static constexpr int NCHILD = 1 << D; //!< 8 for D=3, 4 for D=2.
  std::array<size_t, NCHILD> child; /*!< \brief Indices of child cells. */
  double *M;
  double *L;
  std::vector<size_t> leaf; /*!< \brief Indices of particles in the cell. */
  /*! \brief Start of this cell's particles in the tree's Morton-ordered
      arrays. A leaf's particles occupy [body_offset, body_offset + nleaf).
      Meaningful for leaf cells only. */
  size_t body_offset;
  std::array<double, D> centre; /*!< \brief Coordinates of cell centre. */
  double r; /*!< \brief Radius of cell
                Must be sufficiently large for the root cell to bound the
                particles.

                Note: I may change this in future so it is calculated
                in build_tree rather than user specified.
                */
  double rmax;
  size_t parent; /*!< \brief Index of parent cell of this cell. */
  Cell(std::array<double, D> centre, double r, size_t parent, size_t level, size_t ncrit);
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
    this->centre = other.centre;
    this->r = other.r;
    this->parent = other.parent;
    return *this;
  }
};


template <int D>
class Tree {
public:
  size_t order;
  size_t ncrit;
  double theta;
  std::vector<Particle<D>> particles;
  std::vector<Cell<D>> cells;
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

      A tree subdivided by (comparing each of the D coordinates against the
      cell centre) visits cells in Z-order (Morton order for D=3, its 2D
      analogue for D=2), so emitting each leaf's particles in cell order makes
      every leaf a contiguous range -- no key computation needed.

      Coordinates are held BY VALUE, split into D separate SoA arrays -- one
      per axis, none allocated for axes that don't exist (a 2D tree never
      allocates a z array, matching the "don't store values that are always
      zero" principle used for the D=2 (planar) generated operators). The
      Particle class stores double* into caller memory, so the P2P inner loop
      previously went leaf[p] -> index -> pointer -> scattered load: three
      dependent loads per source particle. SoA gives the vectorised kernel
      unit-stride loads instead of stride-D gathers. */
  std::array<std::vector<double>, D> body;
  //! All-zero, sized to nparticles, used only when D==2: the generated P2P
  //! kernel's signature always takes 3 source coordinate arrays regardless of
  //! D (it is not duplicated for the planar case -- see calculate.cpp), so a
  //! D=2 tree needs *some* array of the right length to hand it as "z". A
  //! single array reused across every call, rather than a per-array zero
  //! shell the way the generated operators would store one, costs O(N) once,
  //! not O(N) per multipole/local coefficient array.
  std::vector<double> body_z_zero;
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

template <int D>
void printTreeParticles(std::vector<Cell<D>> &cells, size_t cell, size_t depth);

template <int D>
void add_child(std::vector<Cell<D>> &cells, size_t octant, size_t p, size_t ncrit);

template <int D>
void split_cell(std::vector<Cell<D>> &cells, std::vector<Particle<D>> &particles, size_t p, size_t ncrit);

//! Build a tree over `nparticles` particles with D-wide positions `pos`
//! (row-major, stride D) and FMMGEN_SOURCESIZE-wide source strengths `mu`.
template <int D>
Tree<D> build_tree(double *pos, double *mu, size_t nparticles, size_t ncrit, size_t order, double theta);
