# distutils: language = c++

from fidimag.atomistic.energy import Energy
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.algorithm cimport sort
cimport numpy as np
import numpy as np
import sys

cdef extern from "utils.hpp":
    size_t Nterms(size_t p)

cdef extern from "operators.h":
    cdef int FMMGEN_MAXORDER

MAXORDER = FMMGEN_MAXORDER

cdef extern from "tree.hpp":
    cdef cppclass Tree:
        size_t order;
        size_t ncrit;
        double theta;
        Tree();
        void compute_field_fmm(double *F)
        void compute_field_exact(double *F)
        void compute_field_bh(double *F)
    Tree build_tree(double *pos, double *mu, size_t nparticles, size_t ncrit, size_t order, double theta)


cdef class FMM:
    cdef public size_t n
    cdef public size_t ncrit
    cdef public double theta
    cdef public size_t order
    cdef public double [:, :] r
    cdef public double [:] mu
    cdef public double [:] mu_s
    cdef public double [:] Mu
    cdef public int calc_type
    cdef Tree tree

    def __cinit__(self, size_t n, size_t ncrit, double theta, size_t order, double [:, :] r, double [:] mu, double [:] mu_s, calc_type=0):
        if order > MAXORDER:
            raise ValueError(f"Order needs to be < {MAXORDER}")
        self.calc_type = calc_type
        # self.particles = vector[Particle]
        self.n = n
        self.ncrit = ncrit
        self.theta = theta
        self.order = order
        # Don't remove these two line, or the memory goes out of scope!
        self.r = r
        self.mu = mu
        self.mu_s = mu_s
        self.Mu = np.zeros(3*self.n)

        self.tree = build_tree(&self.r[0, 0], &self.Mu[0], self.n, self.ncrit, self.order, self.theta)

    cdef _scale(self):
        for i in range(self.n):
            self.Mu[3*i + 0] = self.mu[3*i + 0] * self.mu_s[i]
            self.Mu[3*i + 1] = self.mu[3*i + 1] * self.mu_s[i]
            self.Mu[3*i + 2] = self.mu[3*i + 2] * self.mu_s[i]

    def compute_field(self, double [:] F):
        self._scale()
        if self.calc_type == 0:
            self.tree.compute_field_fmm(&F[0])
        elif self.calc_type == 1:
            self.tree.compute_field_bh(&F[0])

    def compute_field_exact(self, double [:] F):
        self._scale()
        self.tree.compute_field_exact(&F[0])
