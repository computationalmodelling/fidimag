#include "tree.hpp"
#include<cmath>
#include<vector>
#include<array>
#include<algorithm>
#include "utils.hpp"
#include "calculate.hpp"

/*!
* \brief Constructor for Cell class
*
* \param x position of the cell
* \param y position of the cell
* \param z position of the cell
* \param r Radius of the cell
* \param parent Index of parent cell in tree array of cells.
* \param order Multipole expansion order
* \param level Depth of tree at which the cell sits
* \param ncrit Number of particles held before cell splits
*/

Cell::Cell(double x, double y, double z, double r, size_t parent, size_t order, size_t level, size_t ncrit) {
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


// /* Clear expansion array */
// void Cell::clear() {
//   std::fill(M.begin(), M.end(), 0.0);
//   std::fill(L.begin(), L.end(), 0.0);
// }
//
//
// void Cell::resize(size_t order) {
//   this->M.resize(Nterms(order), 0.0);
//   this->L.resize(Nterms(order), 0.0);
// }


/*! \brief Destructor for the Cell class */
Cell::~Cell() {
  #ifdef FMMLIBDEBUG
    std::cout << "Destructor of Cell called" << std::endl;
  #endif
}

/*! \brief Copy constructor for the Cell class */
Cell::Cell(const Cell& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->r = other.r;
    this->rmax = other.rmax;
    this->parent = other.parent;
    this->level = other.level;
    this->child = other.child;
    // std::copy(other.M.begin(), other.M.end(), std::back_inserter(this->M));
    // std::copy(other.L.begin(), other.L.end(), std::back_inserter(this->L));
    std::copy(other.leaf.begin(), other.leaf.end(), std::back_inserter(this->leaf));
    std::copy(other.child.begin(), other.child.end(), std::back_inserter(this->child));
    this->nleaf = other.nleaf;
    this->nchild = other.nchild;
}

/*! \brief Move constructor for the Cell class */
Cell::Cell(Cell&& other) {
  // Move Constructor
  // The omp_destroy_lock function gets called
  // in the destructor. That means we can't naively set other.lock
  // to a null pointer here, but there's no other appropriate
  // place to put the lock destructor call. So in the destructor, before
  // calling omp_destroy_lock(this->lock), first we check if it's a nullptr
  // before we deallocate the memory.
  this->x = other.x;
  this->y = other.y;
  this->z = other.z;
  this->r = other.r;
  this->rmax = other.rmax;
  this->parent = other.parent;
  this->level = other.level;
  this->child = other.child;
  //this->M = std::move(other.M);
  //this->L = std::move(other.L);
  this->M = other.M;
  this->L = other.L;
  this->leaf = other.leaf;
  this->nleaf = other.nleaf;
  this->nchild = other.nchild;

  // other.M.clear();
  // other.L.clear();
  other.leaf.clear();
  other.child.clear();
}

/*! \brief Print function to show particles and their locations in the tree
* \param cells Reference to std::vector containing all cells in a tree.
* \param cell Cell to start from - user should supply '0' as the function is recursive.
* \param depth Depth of the tree - user should supply '0'
*/
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


/*! \brief Given a cell, add a child to it.
* \param cells Reference to std::vector containing all cells in a tree.
* \param octant Octant at which the new cell should be created.
* \param p The parent cell of the new cell.
* \param ncrit The maximum number of particles in a cell before it splits.
* \param order Order of multipole expansions.
*/
void add_child(std::vector<Cell> &cells, int octant, size_t p, size_t ncrit, size_t order) {
    int c = cells.size();
    // Do not change octant to size_t - otherwise the calculation
    // of x, y, z position through bit masking is *not* correct.
    double r = cells[p].r / 2.0;
    double x = cells[p].x + r * ((octant & 1) * 2 - 1);
    double y = cells[p].y + r * ((octant & 2) - 1);
    double z = cells[p].z + r * ((octant & 4) / 2 - 1);
    size_t parent = p;
    size_t level = cells[p].level + 1;
    cells.push_back(Cell(x, y, z, r, parent, order, level, ncrit));
    cells[p].child[octant] = c;
    cells[c].nleaf = 0;
    cells[p].nchild = (cells[p].nchild | (1 << octant));
}

/*! \brief Splits a cell
* When a cell holds more than ncrit particles, the cell must be split.
* Children are added to thec cells list if they have not already been created,
* and particles are reassigned to these child cells.
*
* \param cells Reference to std::vector containing all cells in a tree.
* \param particles Reference to std::vector containing all particles in a tree.
* \param p The cell which is to be split
* \param ncrit The maximum number of particles in a cell before it splits.
* \param order Order of multipole expansions.

*/
void split_cell(std::vector<Cell> &cells, std::vector<Particle> &particles, size_t p, size_t ncrit, size_t order) {
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


/*! \brief Builds a tree from a bounding cell and a set of particles.
* \param particles Reference to a std::vector of Particle objects.
* \param root The root cell of the system.
* \param ncrit The maximum number of particles in a cell before it splits.
* \param order Order of multipole expansions.
* \returns std::vector of Cells that together are the octree.
*/
Tree build_tree(double *pos, double *mu, size_t nparticles, size_t ncrit, size_t order, double theta) {
  // Create particles list for convenience

  std::vector<Particle> particles(nparticles);
  for(size_t i = 0; i < nparticles; i++) {
    particles[i].r = &pos[3*i];
    particles[i].mu = &mu[3*i];
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
  std:: cout << "Building Tree: Avg pos = (" << xavg << ", " << yavg << ", " << zavg << ")" << std::endl;

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
  double r = (xmax > ymax ? (xmax > zmax? xmax: zmax): (ymax > zmax ? ymax: zmax)) * 1.001;
  std::cout << "xmax = " << xmax << ", ymax = " << ymax << ", zmax = " << zmax << ", rmax = " << r << std::endl;
  auto root = Cell(xavg, yavg, zavg, r, 0, order, 0, ncrit);

  cells.push_back(root);
  for(size_t i = 0; i < particles.size(); i++) {
    curr = 0;
    while (cells[curr].nleaf >= ncrit) {
      cells[curr].nleaf += 1;
      octant = (particles[i].r[0] > cells[curr].x) + ((particles[i].r[1] > cells[curr].y) << 1) + ((particles[i].r[2] > cells[curr].z) << 2);
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


  // Now create tree object, and set properties.
  // Choosing a very simple data type here.
  Tree tree;
  tree.theta = theta;
  tree.ncrit = ncrit;
  tree.order = order;
  tree.cells = cells;
  tree.particles = particles;

  // Create interaction lists, and sort M2L list for cache efficiency.
  interact_dehnen_lazy(0, 0, tree.cells, particles, theta, order, ncrit, tree.M2L_list, tree.P2P_list);
  std::sort(tree.M2L_list.begin(), tree.M2L_list.end(),
         [](std::pair<size_t, size_t> &left, std::pair<size_t, size_t> &right) {
              return left.first < right.first;
             }
         );


  // Create memory into which each cell can point for the multipole arrays.
  tree.M.resize(tree.cells.size() * (Nterms(order) - 1), 0.0);
  tree.L.resize(tree.cells.size() * Nterms(order - 1), 0.0);
  for(size_t i = 0; i < tree.cells.size(); i++) {
    tree.cells[i].M = &tree.M[i*(Nterms(order) - Nterms(0))];
    tree.cells[i].L = &tree.L[i*(Nterms(order - 1))];
  }
  return tree;
}

void Tree::clear_expansions() {
  std::fill(M.begin(), M.end(), 0);
  std::fill(L.begin(), L.end(), 0);
}

void Tree::compute_field(double *F) {
  for(size_t i = 0; i < 3*particles.size(); i++) {
    F[i] = 0.0;
  }
  clear_expansions();
  #pragma omp parallel
  evaluate_P2M(particles, cells, 0, ncrit, order);

  evaluate_M2M(particles, cells, order);
  #pragma omp barrier
  #pragma omp parallel
  {
    evaluate_M2L_lazy(cells, M2L_list, order);
    evaluate_P2P_lazy(cells, particles,P2P_list, F);
    #pragma omp barrier
    evaluate_L2L(cells, order);
    #pragma omp barrier
    evaluate_L2P(particles, cells, F, ncrit, order);
  }
}

void Tree::compute_field_exact(double *F) {
  evaluate_direct(particles, F, particles.size());
}