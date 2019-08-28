#ifndef FMMLIB_TREE_HPP
#define FMMLIB_TREE_HPP

#include<iostream>
#include<cmath>
#include<vector>
#include<array>
#include "utils.hpp"

class Particle {
public:
  double x, y, z;
  double mux, muy, muz;
  Particle(double x, double y, double z, double mux, double muy, double muz) : x(x), y(y), z(z), mux(mux), muy(muy), muz(muz) {}
};

struct Cell {
public:
  unsigned int nleaf;
  unsigned int nchild;
  unsigned int level;
  std::array<unsigned int, 8> child;
  std::vector<double> Mx;
  std::vector<double> My;
  std::vector<double> Mz;
  std::vector<double> Lx;
  std::vector<double> Ly;
  std::vector<double> Lz;
  std::vector<unsigned int> leaf;
  double x, y, z, r;
  unsigned int parent;
  Cell(double x, double y, double z, double r, unsigned int parent, unsigned int order, unsigned int level, unsigned int ncrit);
};

void printTreeParticles(std::vector<Cell> &cells, unsigned int cell, unsigned int depth);

void add_child(std::vector<Cell> &cells, int octant, unsigned int p, unsigned int ncrit, unsigned int order);


void split_cell(std::vector<Cell> &cells, std::vector<Particle> &particles, unsigned int p, unsigned int ncrit, unsigned int order);

std::vector<Cell> build_tree(std::vector<Particle> &particles, Cell &root, unsigned int ncrit, unsigned int order);

#endif