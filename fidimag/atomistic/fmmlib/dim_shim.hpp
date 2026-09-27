#pragma once
#include "tree.hpp"

// Cython cannot bind a non-type (int) C++ template parameter directly, so
// this gives it concrete names to bind to instead. Tree<D> is already
// explicitly instantiated for D=2 and D=3 at the bottom of tree.cpp.
typedef Tree<2> Tree2;
typedef Tree<3> Tree3;

inline Tree3 build_tree_3d(double *pos, double *mu, size_t nparticles,
                           size_t ncrit, size_t order, double theta) {
  return build_tree<3>(pos, mu, nparticles, ncrit, order, theta);
}

inline Tree2 build_tree_2d(double *pos, double *mu, size_t nparticles,
                           size_t ncrit, size_t order, double theta) {
  return build_tree<2>(pos, mu, nparticles, ncrit, order, theta);
}
