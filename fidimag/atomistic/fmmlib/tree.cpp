#include<cmath>
#include<vector>
#include<array>
#include "utils.hpp"
#include "tree.hpp"



Particle::Particle(): r(nullptr), mu(nullptr) {}



Cell::Cell(){}




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
    #ifdef FMMLIBDEBUG
      std::cout << "Cell("<<x<<","<<y<<","<<z<<","<<r<<","<<parent<<","<<order<<","<<level<<","<<ncrit<<")"<< std::endl;
    #endif
    this->x = x;
    this->y = y;
    this->z = z;
    this->r = r;
    this->rmax = sqrt(3*r);
    this->parent = parent;
    this->level = level;
    this->child.resize(8, 0);
    // std::cout << "Nterms(order) = " << Nterms(order) << std::endl;
    // std::cout << "Nterms(0) = " << Nterms(0) << std::endl;
    // std::cout << "M.size() = " << Nterms(order) - Nterms(0) << std::endl;
    // this->M.resize(Nterms(order) - Nterms(0), 0.0);
    // this->L.resize(Nterms(order - 1), 0.0);
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
    // of x, y, z position is not correct.
    double r = cells[p].r / 2.0;
    //std::cout << "Adding child to cell p "<< p<< "Parent(" << cells[p].x << "," << cells[p].y << "," << cells[p].z << "," << cells[p].r << ")" << std::endl;

    double x = cells[p].x + r * ((octant & 1) * 2 - 1);
    double y = cells[p].y + r * ((octant & 2) - 1);
    double z = cells[p].z + r * ((octant & 4) / 2 - 1);
    #ifdef FMMLIBDEBUG
        std::cout << "Creating Cell("<<x<<","<<y<<","<<z<<","<<r<<")"<<std::endl;
    #endif
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
  //std::cout << "Printing split_cell" << std::endl;
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
std::vector<Cell> build_tree(std::vector<Particle> &particles, Cell &root, size_t ncrit, size_t order) {
  std::cout << "Building tree" << std::endl;
  std::vector<Cell> cells;
  size_t curr;
  int octant;
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
  std::cout << "Ncells = " << cells.size() << std::endl;
  return cells;
}


/*! \brief Sets multipole expansions to zero.
* \param cells Reference to std::vector containing all cells in a tree.
*/
// void clear_expansions(std::vector<Cell> cells) {
//   for(size_t c = 0; c < cells.size(); c++) {
//     std::fill(cells[c].M.begin(), cells[c].M.end(), 0.0);
//     std::fill(cells[c].L.begin(), cells[c].L.end(), 0.0);
//   }
// }
