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

Cell::Cell(double x, double y, double z, double r, size_t parent, size_t level, size_t ncrit) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;
    this->rmax = sqrt(3*r*r);
    this->parent = parent;
    this->level = level;
    this->child.resize(8, 0);
    this->leaf.resize(ncrit, 0);
    this->nleaf = 0;
    this->nchild = 0;
}

Cell::~Cell() {
  #ifdef FMMLIBDEBUG
    std::cout << "Destructor of Cell called" << std::endl;
  #endif
}

Cell::Cell(const Cell& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->r = other.r;
    this->rmax = other.rmax;
    this->parent = other.parent;
    this->level = other.level;
    this->child = other.child;
    std::copy(other.leaf.begin(), other.leaf.end(), std::back_inserter(this->leaf));
    // Removed duplicate child copy - line 38 already copied it
    this->nleaf = other.nleaf;
    this->nchild = other.nchild;
}

Cell::Cell(Cell&& other) {
  this->x = other.x;
  this->y = other.y;
  this->z = other.z;
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
  other.child.clear();
}

void printTreeParticles(std::vector<Cell> &cells, size_t cell, size_t depth) {
  for(size_t i = 0; i < depth; i++) {
    std::cout << "         ";
  }
  std::cout << cell << " (" << cells[cell].x << ","<< cells[cell].y << "," << cells[cell].z << ") : (";
  size_t nchild = 0;
  for(size_t octant = 0; octant < 8; octant++) {
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
  for(size_t octant = 0; octant < 8; octant++) {
    if (cells[cell].nchild & (1 << octant)) {
      printTreeParticles(cells, cells[cell].child[octant], depth + 1);
    }
  }
}

void add_child(std::vector<Cell> &cells, int octant, size_t p, size_t ncrit) {
    int c = cells.size();
    // Do not change octant to size_t - otherwise the calculation
    // of x, y, z position through bit masking is *not* correct.
    double r = cells[p].r / 2.0;
    double x = cells[p].x + r * ((octant & 1) * 2 - 1);
    double y = cells[p].y + r * ((octant & 2) - 1);
    double z = cells[p].z + r * ((octant & 4) / 2 - 1);
    size_t parent = p;
    size_t level = cells[p].level + 1;
    cells.push_back(Cell(x, y, z, r, parent, level, ncrit));
    cells[p].child[octant] = c;
    cells[c].nleaf = 0;
    cells[p].nchild = (cells[p].nchild | (1 << octant));
}

void split_cell(std::vector<Cell> &cells, std::vector<Particle> &particles, size_t p, size_t ncrit) {
  size_t l, c;
  // Do not change octant to size_t - otherwise the calculation
  // of x, y, z position in add_child is not correct!
  int octant;
  for(size_t i = 0; i < cells[p].leaf.size(); i++) {
    l = cells[p].leaf[i];
    octant = (particles[l].r[0] > cells[p].x) +
      ((particles[l].r[1] > cells[p].y) << 1) +
      ((particles[l].r[2] > cells[p].z) << 2);

    if (!((cells[p].nchild) & (1 << octant))) {
      add_child(cells, octant, p, ncrit);
    }
    c = cells[p].child[octant];
    cells[c].leaf[cells[c].nleaf] = l;
    cells[c].nleaf += 1;
    if (cells[c].nleaf >= ncrit) {
      split_cell(cells, particles, c, ncrit);
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

Tree build_tree(double *pos, double *S, size_t nparticles, size_t ncrit, size_t order, double theta) {
  // Create particles list for convenience

  std::vector<Particle> particles(nparticles);
  for(size_t i = 0; i < nparticles; i++) {
    particles[i].r = &pos[3*i];
    particles[i].S = &S[FMMGEN_SOURCESIZE*i];
  }

  // Now create cells list
  std::vector<Cell> cells;
  size_t curr;
  int octant;

  // Compute average position
  double xavg = 0;
  double yavg = 0;
  double zavg = 0;
  for(size_t i = 0; i < particles.size(); i++) {
    xavg += particles[i].r[0];
    yavg += particles[i].r[1];
    zavg += particles[i].r[2];
  }

  xavg /= particles.size();
  yavg /= particles.size();
  zavg /= particles.size();
  #ifdef FMMLIBDEBUG
  std:: cout << "Building Tree: Avg pos = (" << xavg << ", " << yavg << ", " << zavg << ")" << std::endl;
  #endif
  double xmax = 0;
  double ymax = 0;
  double zmax = 0;

  for(size_t i = 0; i < particles.size(); i++) {
    double x = std::abs(particles[i].r[0] - xavg);
    double y = std::abs(particles[i].r[1] - yavg);
    double z = std::abs(particles[i].r[2] - zavg);

    if (x > xmax)
      xmax = x;
    if (y > ymax)
      ymax = y;
    if (z > zmax)
      zmax = z;
  }

  // if xmax > ymax
  //    then if xmax > zmax, return xmax
  //         else zmax
  // else if ymax > zmax, return ymax
  // else return zmax
  // * 1.001 so that cell slightly bigger than furthest away particle.
  // std::cout << "xmax = " << xmax << ", ymax = " << ymax << ", zmax = " << zmax << ", rmax = " << r << std::endl;

  double r = (xmax > ymax ? (xmax > zmax? xmax: zmax): (ymax > zmax ? ymax: zmax)) * 1.001;
  auto root = Cell(xavg, yavg, zavg, r, 0, 0, ncrit);

  // Reserve capacity to avoid reallocations during tree construction
  // Estimate: ~8 cells per ncrit particles
  cells.reserve(nparticles / ncrit * 8);
  cells.push_back(root);
  for(size_t i = 0; i < particles.size(); i++) {
    curr = 0;
    while (cells[curr].nleaf >= ncrit) {
      cells[curr].nleaf += 1;
      octant = (particles[i].r[0] > cells[curr].x) + ((particles[i].r[1] > cells[curr].y) << 1) + ((particles[i].r[2] > cells[curr].z) << 2);
      if (!(cells[curr].nchild & (1 << octant))) {
        add_child(cells, octant, curr, ncrit);
      }
      curr = cells[curr].child[octant];
    }
    cells[curr].leaf[cells[curr].nleaf] = i;
    cells[curr].nleaf += 1;
    if (cells[curr].nleaf >= ncrit) {
      split_cell(cells, particles, curr, ncrit);
    }
  }


  // Now create tree object, and set properties.
  // Choosing a very simple data type here.
  Tree tree;
  tree.theta = theta;
  tree.ncrit = ncrit;
  tree.order = order;
  tree.cells = cells;
  tree.particles = particles;

  interact_dehnen_lazy(0, 0, tree.cells, particles, theta, order, ncrit, tree.M2L_list, tree.P2P_list);

  // Permute particles into tree order and copy them into flat, by-value
  // arrays. See Tree::body_x for why this matters to the P2P kernel.
  {
    const size_t np = particles.size();
    tree.body_x.resize(np); tree.body_y.resize(np); tree.body_z.resize(np);
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
        tree.body_x[off] = particles[l].r[0];
        tree.body_y[off] = particles[l].r[1];
        tree.body_z[off] = particles[l].r[2];
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
void Tree::clear_M() {
  const size_t n = M.size();
  double *const m = M.data();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; i++) m[i] = 0.0;
}

void Tree::clear_L() {
  const size_t n = L.size();
  double *const l = L.data();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; i++) l[i] = 0.0;
}

//! Zero the persistent Morton-ordered accumulator. See Tree::Fm.
void Tree::clear_Fm() {
  const size_t n = Fm.size();
  double *const f = Fm.data();
  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; i++) f[i] = 0.0;
}

//! See the declaration in tree.hpp for why this runs every solve rather
//! than only once in build_tree.
void Tree::refresh_sources() {
  const size_t n = body_perm.size();
  #pragma omp parallel for schedule(static)
  for (size_t off = 0; off < n; off++) {
    const double *const S = particles[body_perm[off]].S;
    for (int d = 0; d < FMMGEN_SOURCESIZE; d++)
      body_S[FMMGEN_SOURCESIZE*off + d] = S[d];
  }
}

void Tree::compute_field_fmm(double *F) {
  // Computed in Morton order, then scattered once into the caller's ordering.
  // Keeping the whole solve in Morton order means every stage writes
  // contiguously; the only permuted access in the entire method is the single
  // pass at the end.
  refresh_sources();
  clear_Fm();
  clear_M();
  clear_L();
  #pragma omp parallel
  {
    evaluate_P2M(cells, body_x.data(), body_y.data(), body_z.data(),
                 body_S.data(), ncrit, order);
  }
  #pragma omp parallel
  {
    evaluate_M2M(cells, levels, order);
  }
  #ifdef FMMLIBDEBUG
  M_sanity_check(cells);
  #endif
  #pragma omp parallel
  {
    evaluate_M2L_lazy(cells, M2L_list, M2L_group, order);
    evaluate_P2P_lazy(cells, body_x.data(), body_y.data(), body_z.data(),
                      body_S.data(), P2P_list, P2P_group, Fm.data());
  }
  #pragma omp parallel
  {
    evaluate_L2L(cells, levels, order);
  }
  #pragma omp parallel
  {
    evaluate_L2P(cells, body_x.data(), body_y.data(), body_z.data(),
                 Fm.data(), ncrit, order);
  }
  scatter_output(Fm.data(), F);
}

void Tree::compute_field_bh(double *F) {
  refresh_sources();
  clear_Fm();
  clear_M();
  #pragma omp parallel
  {
    evaluate_P2M(cells, body_x.data(), body_y.data(), body_z.data(),
                 body_S.data(), ncrit, order);
  }
  #pragma omp parallel
  {
    evaluate_M2M(cells, levels, order);
  }
  #ifdef FMMLIBDEBUG
  M_sanity_check(cells);
  #endif
  #pragma omp parallel for schedule(runtime)
  for (size_t m = 0; m < particles.size(); m++) {
    evaluate_M2P_and_P2P(body_x.data(), body_y.data(), body_z.data(),
                         body_S.data(), 0, m, cells, Fm.data(), ncrit, theta, order);
  }
  scatter_output(Fm.data(), F);
}

void Tree::compute_field_exact(double *F) {
  // Keeps its own scratch buffer rather than reusing the persistent Tree::Fm:
  // this is the reference solution the approximate methods are measured
  // against, so it must not share state with them.
  refresh_sources();
  std::vector<double> Fdirect(FMMGEN_OUTPUTSIZE * particles.size(), 0.0);
  evaluate_direct(body_x.data(), body_y.data(), body_z.data(),
                  body_S.data(), Fdirect.data(), particles.size());
  scatter_output(Fdirect.data(), F);
}

void Tree::scatter_output(const double *Fm, double *F) {
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
