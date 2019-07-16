# distutils: language = c++

from fidimag.atomistic.energy import Energy
from libcpp.vector cimport vector
from libcpp.utility cimport pair

cimport numpy as np
import numpy as np

cdef extern from "utils.hpp":
    size_t Nterms(size_t p)

cdef extern from "tree.hpp":
    cdef cppclass Particle:
        double *r
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

    void evaluate_direct(vector[Particle] particles, 
                         double *F,
                         size_t n)





cdef class FMM:
    cdef size_t n, ncrit, order
    cdef vector[Particle] particles
    cdef vector[Cell] cells
    cdef size_t i
    cdef Particle p
    cdef Cell root
    cdef vector[double] M
    cdef vector[double] L
    cdef size_t Msize
    cdef size_t Lsize
    cdef public double [:, :] coords

    def __cinit__(self, size_t n, size_t ncrit, size_t order, double [:, :] coords, double [:] mu):
        # self.particles = vector[Particle]

        # Don't remove this line, or the memory goes out of scope!
        self.coords = coords
        self.ncrit = ncrit
        self.order = order
        print('FMM Order = {}'.format(order))
        xs = np.asarray(self.coords[:, 0])
        ys = np.asarray(self.coords[:, 1])
        zs = np.asarray(self.coords[:, 2])

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
            self.p = Particle(&self.coords[i, 0],
                              &mu[3*i])
            self.particles.push_back(self.p)
        self.cells = build_tree(self.particles, self.root, ncrit, order)
        print(f"Cython: cells.size() {self.cells.size()}")


        print(f"r[0] = {self.particles[0].r[0]}, {self.particles[0].r[1]}, {self.particles[0].r[2]}")

        print(f"mu[0] = {self.particles[0].mu[0]}, {self.particles[0].mu[1]}, {self.particles[0].mu[2]}")

        # Need to allocate memory for the M and L arrays.
        #vector[double] M(cells.size() * (Nterms(order) - Nterms(0)), 0.0)
        #vector[double] L(cells.size() * Nterms(order - 1), 0.0)

        self.Msize = Nterms(order) - Nterms(0)
        self.Lsize = Nterms(order - 1)
        self.M = vector[double](self.cells.size() * self.Msize, 0.0)
        self.L = vector[double](self.cells.size() * self.Lsize, 0.0)

        for i in range(self.cells.size()):
            self.cells[i].M = &self.M[i*self.Msize]
            self.cells[i].L = &self.L[i*self.Lsize]


    def P2M(self):
        evaluate_P2M(self.particles, self.cells, 0, self.ncrit, self.order)

    def M2M(self):
        evaluate_M2M(self.particles, self.cells, self.order)


    def compute_field(self, double theta, double [:] F):
        F[:] = 0
        print('compute field')
        print("P2M starting")
        evaluate_P2M(self.particles, self.cells, 0, self.ncrit, self.order)
        print("M2M starting")
        evaluate_M2M(self.particles, self.cells, self.order)

        print(f"mu[0] = {self.particles[0].mu[0]}, {self.particles[0].mu[1]}, {self.particles[0].mu[2]}")

        print(np.max(self.M))

        print("interact_dehnen starting")
        interact_dehnen(0, 0, self.cells, self.particles, theta, self.order, self.ncrit, &F[0])
        print("L2L starting")
        evaluate_L2L(self.cells, self.order)
        print("L2P starting")
        evaluate_L2P(self.particles, self.cells, &F[0], self.ncrit, self.order)
        print(F[0])

    def compute_field_exact(self, double [:] F_exact):
        evaluate_direct(self.particles, &F_exact[0], self.n)

    def numcells(self):
        return self.cells.size()

    def printM(self, cell, upto=3):
        if cell > self.cells.size():
            raise ValueError(f"Cell index must be less than {self.cells.size()}")
        else:
            for j in range(upto):
                print f"M[{j}] = {self.cells[cell].M[j]}"

    def printL(self, cell):
        if cell > self.cells.size():
            raise ValueError(f"Cell index must be less than {self.cells.size()}")
        else:
            for j in range(self.Lsize):
                print f"L[{j}] = {self.L[cell*self.Lsize + j]}"

    def part(self, i):
        return (np.array([self.particles[i].r[0], self.particles[i].r[1], self.particles[2].r[2]]),
                np.array([self.particles[i].mu[0], self.particles[i].mu[1], self.particles[2].mu[2]]))
