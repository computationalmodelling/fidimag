#include "tree.hpp"
#include "variant.hpp"
#include <stdexcept>
#include <string>
#include<cmath>
#include<vector>
#include<array>
#include<algorithm>
#include "utils.hpp"
#include "calculate.hpp"
#include "operators.h"

template <int D>
Cell<D>::Cell(std::array<double, D> centre, double r, size_t parent, size_t level, size_t ncrit) {
    this->centre = centre;
    this->r = r;
    this->rmax = sqrt(D * r * r);
    this->parent = parent;
    this->level = level;
    this->child.fill(0);
    this->leaf.resize(ncrit, 0);
    this->nleaf = 0;
    this->nchild = 0;
}

template <int D>
Cell<D>::~Cell() {
  #ifdef FMMLIBDEBUG
    std::cout << "Destructor of Cell called" << std::endl;
  #endif
}

template <int D>
Cell<D>::Cell(const Cell& other) {
    this->centre = other.centre;
    this->r = other.r;
    this->rmax = other.rmax;
    this->parent = other.parent;
    this->level = other.level;
    this->child = other.child;
    std::copy(other.leaf.begin(), other.leaf.end(), std::back_inserter(this->leaf));
    this->nleaf = other.nleaf;
    this->nchild = other.nchild;
}

template <int D>
Cell<D>::Cell(Cell&& other) {
  this->centre = other.centre;
  this->r = other.r;
  this->rmax = other.rmax;
  this->parent = other.parent;
  this->level = other.level;
  this->child = other.child;
  this->M = other.M;
  this->L = other.L;
  this->leaf = other.leaf;
  this->nleaf = other.nleaf;
  this->nchild = other.nchild;
  other.leaf.clear();
}

template <int D>
void printTreeParticles(std::vector<Cell<D>> &cells, size_t cell, size_t depth) {
  for(size_t i = 0; i < depth; i++) {
    std::cout << "         ";
  }
  std::cout << cell << " (";
  for (int d = 0; d < D; d++) {
    std::cout << cells[cell].centre[d];
    if (d + 1 < D) std::cout << ",";
  }
  std::cout << ") : (";
  size_t nchild = 0;
  for(size_t octant = 0; octant < Cell<D>::NCHILD; octant++) {
    if (cells[cell].nchild & (1 << octant)) {
      nchild += 1;
    }
  }

  if (nchild == 0) {
    for(size_t i = 0; i < cells[cell].nleaf; i++) {
      std::cout << cells[cell].leaf[i];
      if (i != (cells[cell].nleaf - 1)) {
        std::cout << ",";
      }
    }
  }
  std::cout << ")" << std::endl;
  for(size_t octant = 0; octant < Cell<D>::NCHILD; octant++) {
    if (cells[cell].nchild & (1 << octant)) {
      printTreeParticles(cells, cells[cell].child[octant], depth + 1);
    }
  }
}

// Which of the 2^D children a point falls into, relative to a cell centre:
// bit d is set if the point's coordinate d is greater than the cell centre's.
template <int D>
static inline size_t octant_of(const double *r, const std::array<double, D> &centre) {
  size_t octant = 0;
  for (int d = 0; d < D; d++) {
    octant |= (size_t)(r[d] > centre[d]) << d;
  }
  return octant;
}

template <int D>
void add_child(std::vector<Cell<D>> &cells, size_t octant, size_t p, size_t ncrit) {
    size_t c = cells.size();
    double r = cells[p].r / 2.0;
    std::array<double, D> centre;
    for (int d = 0; d < D; d++) {
      // The bit extract must be forced to signed arithmetic: octant is
      // unsigned, and "0 * 2 - 1" computed in an unsigned type underflows to
      // a huge positive value instead of -1. See split_cell's comment on the
      // same trap in the original (per-axis) formula this generalises.
      int bit = (int)((octant >> d) & 1);
      centre[d] = cells[p].centre[d] + r * (bit * 2 - 1);
    }
    size_t parent = p;
    size_t level = cells[p].level + 1;
    cells.push_back(Cell<D>(centre, r, parent, level, ncrit));
    cells[p].child[octant] = c;
    cells[c].nleaf = 0;
    cells[p].nchild = (cells[p].nchild | (1 << octant));
}

template <int D>
void split_cell(std::vector<Cell<D>> &cells, std::vector<Particle<D>> &particles, size_t p, size_t ncrit) {
  size_t l, c;
  size_t octant;
  for(size_t i = 0; i < cells[p].leaf.size(); i++) {
    l = cells[p].leaf[i];
    octant = octant_of<D>(particles[l].r, cells[p].centre);

    if (!((cells[p].nchild) & (1 << octant))) {
      add_child<D>(cells, octant, p, ncrit);
    }
    c = cells[p].child[octant];
    cells[c].leaf[cells[c].nleaf] = l;
    cells[c].nleaf += 1;
    if (cells[c].nleaf >= ncrit) {
      split_cell<D>(cells, particles, c, ncrit);
  }
  }
}

// Stable counting sort of an interaction list by target cell, also returning
// the bucket offsets: entries for target A end up in [group[A], group[A+1]).
static void group_by_target(std::vector<Interaction> &list,
                            std::vector<size_t> &group, size_t ncells) {
  const size_t n = list.size();
  group.assign(ncells + 1, 0);
  for (size_t i = 0; i < n; i++) group[list[i].target + 1]++;
  for (size_t c = 0; c < ncells; c++) group[c + 1] += group[c];

  std::vector<size_t> cursor(group.begin(), group.end() - 1);
  std::vector<Interaction> out(n);
  for (size_t i = 0; i < n; i++) out[cursor[list[i].target]++] = list[i];
  list.swap(out);
}

static void check_grouped(const std::vector<Interaction> &list, const char *name) {
  for (size_t i = 1; i < list.size(); i++) {
    if (list[i].target < list[i-1].target) {
      throw std::runtime_error(std::string(name) + " not grouped by target: "
                               "grouped evaluation would race on the target accumulator");
    }
  }
}

template <int D>
Tree<D> build_tree(double *pos, double *S, size_t nparticles, size_t ncrit, size_t order, double theta) {
  // Create particles list for convenience

  std::vector<Particle<D>> particles(nparticles);
  for(size_t i = 0; i < nparticles; i++) {
    particles[i].r = &pos[D*i];
    particles[i].S = &S[FMMGEN_SOURCESIZE*i];
  }

  // Now create cells list
  std::vector<Cell<D>> cells;
  size_t curr, octant;

  // Compute average position
  std::array<double, D> avg{};
  for(size_t i = 0; i < particles.size(); i++) {
    for (int d = 0; d < D; d++) avg[d] += particles[i].r[d];
  }
  for (int d = 0; d < D; d++) avg[d] /= particles.size();
  #ifdef FMMLIBDEBUG
  std::cout << "Building Tree: Avg pos = (";
  for (int d = 0; d < D; d++) std::cout << avg[d] << (d+1<D ? ", " : "");
  std::cout << ")" << std::endl;
  #endif

  std::array<double, D> maxext{};
  for(size_t i = 0; i < particles.size(); i++) {
    for (int d = 0; d < D; d++) {
      double v = std::abs(particles[i].r[d] - avg[d]);
      if (v > maxext[d]) maxext[d] = v;
    }
  }

  // * 1.001 so that cell slightly bigger than furthest away particle.
  double r = 0.0;
  for (int d = 0; d < D; d++) r = std::max(r, maxext[d]);
  r *= 1.001;
  auto root = Cell<D>(avg, r, 0, 0, ncrit);

  // Reserve capacity to avoid reallocations during tree construction
  // Estimate: ~(2^D) cells per ncrit particles
  cells.reserve(nparticles / ncrit * Cell<D>::NCHILD);
  cells.push_back(root);
  for(size_t i = 0; i < particles.size(); i++) {
    curr = 0;
    while (cells[curr].nleaf >= ncrit) {
      cells[curr].nleaf += 1;
      octant = octant_of<D>(particles[i].r, cells[curr].centre);
      if (!(cells[curr].nchild & (1 << octant))) {
        add_child<D>(cells, octant, curr, ncrit);
      }
      curr = cells[curr].child[octant];
    }
    cells[curr].leaf[cells[curr].nleaf] = i;
    cells[curr].nleaf += 1;
    if (cells[curr].nleaf >= ncrit) {
      split_cell<D>(cells, particles, curr, ncrit);
    }
  }


  // Now create tree object, and set properties.
  Tree<D> tree;
  tree.theta = theta;
  tree.ncrit = ncrit;
  tree.order = order;
  tree.cells = cells;
  tree.particles = particles;

  interact_dehnen_lazy<D>(0, 0, tree.cells, particles, theta, order, ncrit, tree.M2L_list, tree.P2P_list);

  // Permute particles into tree order and copy them into flat, by-value
  // arrays. See Tree::body for why this matters to the P2P kernel.
  {
    const size_t np = particles.size();
    for (int d = 0; d < D; d++) tree.body[d].resize(np);
    if constexpr (D == 2) tree.body_z_zero.assign(np, 0.0);
    tree.body_S.resize(FMMGEN_SOURCESIZE * np);
    tree.body_perm.resize(np);
    tree.Fm.resize(FMMGEN_OUTPUTSIZE * np);   // persistent; re-zeroed per solve

    size_t off = 0;
    for (size_t c = 0; c < tree.cells.size(); c++) {
      if (tree.cells[c].nleaf >= ncrit) continue;   // internal cell
      tree.cells[c].body_offset = off;
      for (size_t i = 0; i < tree.cells[c].nleaf; i++) {
        const size_t l = tree.cells[c].leaf[i];
        tree.body_perm[off] = l;
        for (int d = 0; d < D; d++) tree.body[d][off] = particles[l].r[d];
        off++;
      }
    }
    if (off != np) {
      throw std::runtime_error("particle permutation did not cover every particle");
    }
    tree.refresh_sources();

    // Cell::leaf has no readers left: every evaluation stage now addresses
    // particles as the contiguous range [body_offset, body_offset + nleaf) in
    // the arrays above. Release it rather than carry a second, redundant
    // representation of the same information (ncrit size_t per cell, and it
    // was allocated for internal cells too).
    for (size_t c = 0; c < tree.cells.size(); c++) {
      std::vector<size_t>().swap(tree.cells[c].leaf);
    }
  }

  // Group both interaction lists by target with a counting sort.
  //
  // The key is a cell index in [0, ncells) -- dense and small (n = 1.8M against
  // k = 4154 cells for a 20k-particle run), so a comparison sort pays ~log2(n)
  // = 21 levels of compare-and-branch for a key it could bucket directly. Two
  // linear passes instead: measured 5.3x faster at 20k and 5.6x at 100k.
  //
  // The exclusive prefix sum IS the group-offset array, so no separate pass is
  // needed to find group boundaries. The sort is stable, so source order within
  // a target follows traversal order rather than depending on tie-breaking,
  // which combined with the atomic-free evaluation makes results reproducible.
  group_by_target(tree.M2L_list, tree.M2L_group, tree.cells.size());
  group_by_target(tree.P2P_list, tree.P2P_group, tree.cells.size());

  // evaluate_M2L_lazy and evaluate_P2P_lazy each give a thread exclusive
  // ownership of one target's accumulator and drop all atomics on that basis.
  // That is only sound if every entry for a target is contiguous. Check it
  // rather than trust it -- getting this wrong is a race, not a clean failure.
  check_grouped(tree.M2L_list, "M2L_list");
  check_grouped(tree.P2P_list, "P2P_list");


  // Create memory into which each cell can point for the multipole arrays.
  // Strides come from the selected variant: the compressed operators use
  // (order+1)^2 coefficients, not Nterms(order). Baked into the cell pointers
  // here, so nothing downstream needs to know which variant is running.
  const size_t ms = fmm->msize(order), ls = fmm->lsize(order);
  tree.M.resize(tree.cells.size() * ms, 0.0);
  tree.L.resize(tree.cells.size() * ls, 0.0);
  for(size_t i = 0; i < tree.cells.size(); i++) {
    tree.cells[i].M = &tree.M[i*ms];
    tree.cells[i].L = &tree.L[i*ls];
  }

  // Precompute tree levels for efficient parallel M2M and L2L
  size_t max_level = 0;
  for(size_t c = 0; c < tree.cells.size(); c++) {
    if (tree.cells[c].level > max_level) {
      max_level = tree.cells[c].level;
    }
  }
  tree.levels.resize(max_level + 1);
  for(size_t c = 0; c < tree.cells.size(); c++) {
    tree.levels[tree.cells[c].level].push_back(c);
  }

  return tree;
}

// M and L are each ncells * Msize(order) doubles -- 16.8 MB apiece at N=1M,
// p=5 -- and were filled by one thread while every other thread waited.
// schedule(static) so each thread re-touches the same pages on every solve.
template <int D>
void Tree<D>::clear_M() {
  const size_t n = M.size();
  double *const m = M.data();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; i++) m[i] = 0.0;
}

template <int D>
void Tree<D>::clear_L() {
  const size_t n = L.size();
  double *const l = L.data();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; i++) l[i] = 0.0;
}

//! Zero the persistent Morton-ordered accumulator. See Tree::Fm.
template <int D>
void Tree<D>::clear_Fm() {
  const size_t n = Fm.size();
  double *const f = Fm.data();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; i++) f[i] = 0.0;
}

//! See the declaration in tree.hpp for why this runs every solve rather
//! than only once in build_tree.
template <int D>
void Tree<D>::refresh_sources() {
  const size_t n = body_perm.size();
  #pragma omp parallel for schedule(static)
  for (size_t off = 0; off < n; off++) {
    const double *const S = particles[body_perm[off]].S;
    for (int d = 0; d < FMMGEN_SOURCESIZE; d++)
      body_S[FMMGEN_SOURCESIZE*off + d] = S[d];
  }
}

// The generated P2P kernel always takes 3 source-coordinate arrays: bz is
// body[2].data() for a real 3D tree, or the persistent zero array for a
// planar (D=2) one. See Tree::body_z_zero.
template <int D>
static inline const double* bz_ptr(const std::array<std::vector<double>, D> &body,
                                    const std::vector<double> &zero) {
  if constexpr (D == 3) return body[2].data();
  else return zero.data();
}

template <int D>
void Tree<D>::compute_field_fmm(double *F) {
  // Computed in Morton order, then scattered once into the caller's ordering.
  // Keeping the whole solve in Morton order means every stage writes
  // contiguously; the only permuted access in the entire method is the single
  // pass at the end.
  refresh_sources();
  clear_Fm();
  clear_M();
  clear_L();
  const double *bx = body[0].data(), *by = body[1].data(), *bz = bz_ptr<D>(body, body_z_zero);
  #pragma omp parallel
  {
    evaluate_P2M<D>(cells, bx, by, bz, body_S.data(), ncrit, order);
  }
  #pragma omp parallel
  {
    evaluate_M2M<D>(cells, levels, order);
  }
  #ifdef FMMLIBDEBUG
  M_sanity_check<D>(cells);
  #endif
  #pragma omp parallel
  {
    evaluate_M2L_lazy<D>(cells, M2L_list, M2L_group, order);
    evaluate_P2P_lazy<D>(cells, bx, by, bz, body_S.data(), P2P_list, P2P_group, Fm.data());
  }
  #pragma omp parallel
  {
    evaluate_L2L<D>(cells, levels, order);
  }
  #pragma omp parallel
  {
    evaluate_L2P<D>(cells, bx, by, bz, Fm.data(), ncrit, order);
  }
  scatter_output(Fm.data(), F);
}

template <int D>
void Tree<D>::compute_field_bh(double *F) {
  refresh_sources();
  clear_Fm();
  clear_M();
  const double *bx = body[0].data(), *by = body[1].data(), *bz = bz_ptr<D>(body, body_z_zero);
  #pragma omp parallel
  {
    evaluate_P2M<D>(cells, bx, by, bz, body_S.data(), ncrit, order);
  }
  #pragma omp parallel
  {
    evaluate_M2M<D>(cells, levels, order);
  }
  #ifdef FMMLIBDEBUG
  M_sanity_check<D>(cells);
  #endif
  #pragma omp parallel for schedule(runtime)
  for (size_t m = 0; m < particles.size(); m++) {
    evaluate_M2P_and_P2P<D>(bx, by, bz, body_S.data(), 0, m, cells, Fm.data(), ncrit, theta, order);
  }
  scatter_output(Fm.data(), F);
}

template <int D>
void Tree<D>::compute_field_exact(double *F) {
  // Keeps its own scratch buffer rather than reusing the persistent Tree::Fm:
  // this is the reference solution the approximate methods are measured
  // against, so it must not share state with them.
  refresh_sources();
  const double *bx = body[0].data(), *by = body[1].data(), *bz = bz_ptr<D>(body, body_z_zero);
  std::vector<double> Fdirect(FMMGEN_OUTPUTSIZE * particles.size(), 0.0);
  evaluate_direct<D>(bx, by, bz, body_S.data(), Fdirect.data(), particles.size());
  scatter_output(Fdirect.data(), F);
}

template <int D>
void Tree<D>::scatter_output(const double *Fm, double *F) {
  // body_perm is a permutation, so every destination index is distinct: this is
  // a pure copy with no accumulation, hence race-free in parallel and unable to
  // change a single digit of the result however it is scheduled.
  //
  // It was the largest serial stage left in the solve: 21.5 ms of 160 ms at
  // N=1M on 96 cores, 13% of the method, and flat in thread count.
  const size_t np = particles.size();
  const size_t *const perm = body_perm.data();
  #pragma omp parallel for schedule(static)
  for (size_t m = 0; m < np; m++) {
    const size_t o = perm[m];
    for (int k = 0; k < FMMGEN_OUTPUTSIZE; k++) {
      F[FMMGEN_OUTPUTSIZE*o + k] = Fm[FMMGEN_OUTPUTSIZE*m + k];
    }
  }
}

// D = 2 (quadtree, planar) and D = 3 (octree) are the only supported
// instantiations -- see the static_assert in Particle<D>.
template class Cell<2>;
template class Cell<3>;
template class Tree<2>;
template class Tree<3>;
template void printTreeParticles<2>(std::vector<Cell<2>> &, size_t, size_t);
template void printTreeParticles<3>(std::vector<Cell<3>> &, size_t, size_t);
template void add_child<2>(std::vector<Cell<2>> &, size_t, size_t, size_t);
template void add_child<3>(std::vector<Cell<3>> &, size_t, size_t, size_t);
template void split_cell<2>(std::vector<Cell<2>> &, std::vector<Particle<2>> &, size_t, size_t);
template void split_cell<3>(std::vector<Cell<3>> &, std::vector<Particle<3>> &, size_t, size_t);
template Tree<2> build_tree<2>(double *, double *, size_t, size_t, size_t, double);
template Tree<3> build_tree<3>(double *, double *, size_t, size_t, size_t, double);
