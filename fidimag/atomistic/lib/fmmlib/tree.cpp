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


Cell::Cell(double x, double y, double z, double r, unsigned int parent, unsigned int order, unsigned int level, unsigned int ncrit) {
    //std::cout << "Cell("<<x<<","<<y<<","<<z<<","<<r<<","<<parent<<","<<order<<","<<level<<","<<ncrit<<")"<< std::endl;
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;
    this->parent = parent;
    this->level = level;
    this->Mx.resize(Nterms(order));
    this->My.resize(Nterms(order));
    this->Mz.resize(Nterms(order));
    for(unsigned int i = 0; i < Mx.size(); i++) {
      this->Mx[i] = 0;
      this->My[i] = 0;
      this->Mz[i] = 0;
    }
    for(unsigned int i = 0; i < leaf.size(); i++) {
      this->leaf[i] = 0;
    }
    this->leaf.resize(ncrit);
    this->nleaf = 0;
    this->nchild = 0;
}


void printTreeParticles(std::vector<Cell> &cells, unsigned int cell, unsigned int depth) {
  if (cell == 0) {
    //std::cout << "Printing particle locations!" << std::endl;
  }
  for(unsigned int i = 0; i < depth; i++) {
    std::cout << "         ";
  }
  std::cout << cell << " (" << cells[cell].x << ","<< cells[cell].y << "," << cells[cell].z << ") : (";
  unsigned int nchild = 0;
  for(unsigned int octant = 0; octant < 8; octant++) {
    if (cells[cell].nchild & (1 << octant)) {
      nchild += 1;
    }
  }

  if (nchild == 0) {
    for(unsigned int i = 0; i < cells[cell].nleaf; i++) {
      std::cout << cells[cell].leaf[i];
      if (i != (cells[cell].nleaf - 1)) {
        std::cout << ",";
      }
    }
  }
  std::cout << ")" << std::endl;
  for(unsigned int octant = 0; octant < 8; octant++) {
    if (cells[cell].nchild & (1 << octant)) {
      printTreeParticles(cells, cells[cell].child[octant], depth + 1);
    }
  }
}

void add_child(std::vector<Cell> &cells, int octant, unsigned int p, unsigned int ncrit, unsigned int order) {
    int c = cells.size();
    // Do not change octant to unsigned int - otherwise the calculation
    // of x, y, z position is not correct.
    double r = cells[p].r / 2.0;
    //std::cout << "Adding child to cell p "<< p<< "Parent(" << cells[p].x << "," << cells[p].y << "," << cells[p].z << "," << cells[p].r << ")" << std::endl;

    double x = cells[p].x + r * ((octant & 1) * 2 - 1);
    double y = cells[p].y + r * ((octant & 2) - 1);
    double z = cells[p].z + r * ((octant & 4) / 2 - 1);
    //std::cout << "Creating Cell("<<x<<","<<y<<","<<z<<","<<r<<")"<<std::endl;
    unsigned int parent = p;
    unsigned int level = cells[p].level + 1;
    cells.push_back(Cell(x, y, z, r, parent, order, level, ncrit));
    cells[p].child[octant] = c;
    cells[c].nleaf = 0;
    cells[p].nchild = (cells[p].nchild | (1 << octant));
}

void split_cell(std::vector<Cell> &cells, std::vector<Particle> &particles, unsigned int p, unsigned int ncrit, unsigned int order) {
  //std::cout << "Printing split_cell" << std::endl;
  unsigned int l, c;
  // Do not change octant to unsigned int - otherwise the calculation
  // of x, y, z position in add_child is not correct!
  int octant;
  for(unsigned int i = 0; i < cells[p].leaf.size(); i++) {
    l = cells[p].leaf[i];
    octant = (particles[l].x > cells[p].x) +
      ((particles[l].y > cells[p].y) << 1) +
      ((particles[l].z > cells[p].z) << 2);

    if (!((cells[p].nchild) & (1 << octant))) {
      add_child(cells, octant, p, ncrit, order);
    }
    c = cells[p].child[octant];
    cells[c].leaf[cells[c].nleaf] = l;
    cells[c].nleaf += 1;
    if (cells[c].nleaf >= ncrit) {
      split_cell(cells, particles, c, ncrit, order);
    }
  }
}

std::vector<Cell> build_tree(std::vector<Particle> &particles, Cell &root, unsigned int ncrit, unsigned int order) {
  std::vector<Cell> cells;
  unsigned int curr;
  int octant;
  cells.push_back(root);
  for(unsigned int i = 0; i < particles.size(); i++) {
    //std::cout << i << std::endl;
    curr = 0;
    while (cells[curr].nleaf >= ncrit) {
      cells[curr].nleaf += 1;
      octant = (particles[i].x > cells[curr].x) + ((particles[i].y > cells[curr].y) << 1) + ((particles[i].z > cells[curr].z) << 2);
      if (!(cells[curr].nchild & (1 << octant))) {
        add_child(cells, octant, curr, ncrit, order);
      }
      curr = cells[curr].child[octant];
    }
    cells[curr].leaf[cells[curr].nleaf] = i;
    cells[curr].nleaf += 1;
    if (cells[curr].nleaf >= ncrit) {
      split_cell(cells, particles, curr, ncrit, order);
    }
  }
  return cells;
}


#endif
