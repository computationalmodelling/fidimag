# distutils: language = c++

from fidimag.atomistic.energy import Energy
from libcpp.vector cimport vector

cimport numpy as np
import numpy as np

cdef extern from "tree.hpp":
    cdef cppclass Particle:
        double x, y, z
        double mux, muy, muz
        Particle()
        Particle(double x, double y, double z, double mux, double muy, double muz)


    cdef cppclass Cell:
        double x, y, z, r, rmax
        vector[double] M
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
    

    void evaluate_approx(vector[Particle] particles,
                         vector[Cell] cells,
                         double *F,
                         unsigned int n_crit,
                         double theta,
                         unsigned int exp_order)

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
            self.p = Particle(coords[i, 0], coords[i, 1], coords[i, 2],
                              mu[3*i+0], mu[3*i+1], mu[3*i+2])
            self.particles.push_back(self.p)

        self.cells = build_tree(self.particles, self.root, ncrit, order)

    def set(self, mu):
        clear_expansions(self.cells)
        for self.i in range(self.n):
            self.particles[self.i].mux = mu[3*self.i+0]
            self.particles[self.i].muy = mu[3*self.i+1]
            self.particles[self.i].muz = mu[3*self.i+2]
        
    def compute_field(self, double theta, double [:] F):
        evaluate_P2M(self.particles, self.cells, 0, self.ncrit, self.order)
        evaluate_M2M(self.particles, self.cells, self.order)
        evaluate_approx(self.particles, self.cells, &F[0], self.ncrit, self.theta, self.order)

    def numcells(self):
        return self.cells.size()


    def getMultipoles(self, i):
        if i >= self.numcells():
            raise ValueError("Cell does not exist")

        else:
            a = np.zeros(self.cells[i].M.size())
            for j in range(len(a)):
                a[j] = self.cells[i].M[j]
            return a
