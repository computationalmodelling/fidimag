#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include<array>
#include "utils.hpp"


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
  Cell(double x, double y, double z, double r, size_t parent, size_t order, size_t level, size_t ncrit);
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
  std::vector<std::pair<size_t, size_t>> M2L_list;
  std::vector<std::pair<size_t, size_t>> P2P_list;
  void compute_field_fmm(double *F);
  void compute_field_bh(double *F);
  void compute_field_exact(double *F);
private:
  void clear_M();
  void clear_L();
};

void printTreeParticles(std::vector<Cell> &cells, size_t cell, size_t depth);

void add_child(std::vector<Cell> &cells, int octant, size_t p, size_t ncrit, size_t order);

void split_cell(std::vector<Cell> &cells, std::vector<Particle> &particles, size_t p, size_t ncrit, size_t order);

Tree build_tree(double *pos, double *mu, size_t nparticles, size_t ncrit, size_t order, double theta);
