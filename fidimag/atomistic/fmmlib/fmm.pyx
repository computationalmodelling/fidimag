# distutils: language = c++

from fidimag.atomistic.energy import Energy
from libcpp.vector cimport vector

cimport numpy as np
import numpy as np

cdef extern from "tree.hpp":
    cdef cppclass Particle:
        double x, y, z
        double *mu
        Particle()
        Particle(double *r, double *mu)


    cdef cppclass Cell:
        double x, y, z, r, rmax
        double *M
        double *L
        Cell(double x, double y, double z, double r, size_t parent, size_t order, size_t level, size_t ncrit)
        Cell()

    vector[Cell] build_tree(vector[Particle] particles, Cell root, size_t ncrit, size_t order)
    void clear_expansions(vector[Cell] cells)



cdef extern from "calculate.hpp":
    void evaluate_P2M(vector[Particle] particles,
                      vector[Cell] cells,
                      size_t cell,
                      size_t ncrit,
                      size_t order)

    void evaluate_M2M(vector[Particle] particles,
                      vector[Cell] cells,
                      size_t order)

    void interact_dehnen(size_t A, size_t B,
                         vector[Cell] cells,
                         vector[Particle] particles,
                         double theta,
                         size_t order,
                         size_t ncrit,
                         double *F)

    void evaluate_L2L(vector[Cell] cells,
                      size_t order)

    void evaluate_L2P(vector[Particle] particles,
                      vector[Cell] cells,
                      double *F,
                      size_t ncrit,
                      size_t order)


cdef class FMM:
    cdef size_t n, ncrit, order
    cdef vector[Particle] particles
    cdef vector[Cell] cells
    cdef size_t i
    cdef Particle p
    cdef Cell root

    def __cinit__(self, size_t n, size_t ncrit, size_t order, double [:, :] coords, double [:] mu):
         # self.particles = vector[Particle]
        self.ncrit = ncrit
        self.order = order
        print('FMM Order = {}'.format(order))
        xs = np.asarray(coords[:, 0])
        ys = np.asarray(coords[:, 1])
        zs = np.asarray(coords[:, 2])

        cdef double rootx = np.mean(xs)
        cdef double rooty = np.mean(ys)
        cdef double rootz = np.mean(zs)

        maxs = np.array([np.abs(xs - rootx),
                        np.abs(ys - rooty),
                        np.abs(zs - rootz)])
        cdef double rmax = np.max(maxs)

        self.root = Cell(rootx, rooty, rootz, rmax, 0, order, 0, ncrit)
        self.n = n

        for i in range(self.n):
            self.p = Particle(&coords[i, 0],
                              &mu[3*i])
            self.particles.push_back(self.p)
        self.cells = build_tree(self.particles, self.root, ncrit, order)

    def compute_field(self, double theta, double [:] F):
        evaluate_P2M(self.particles, self.cells, 0, self.ncrit, self.order)
        evaluate_M2M(self.particles, self.cells, self.order)
        interact_dehnen(0, 0, self.cells, self.particles, theta, self.order, self.ncrit, &F[0])
        evaluate_L2L(self.cells, self.order)
        evaluate_L2P(self.particles, self.cells, &F[0], self.ncrit, self.order)

    # def set(self, m):




    def numcells(self):
        return self.cells.size()
