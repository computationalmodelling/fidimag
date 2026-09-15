# distutils: language = c++

from fidimag.atomistic.energy import Energy
from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.algorithm cimport sort
from libcpp cimport bool
cimport numpy as np
import numpy as np
import sys

cdef extern from "utils.hpp":
    size_t Nterms(size_t p)

cdef extern from "operators.h":
    cdef int FMMGEN_MINORDER
    cdef int FMMGEN_MAXORDER

MINORDER = FMMGEN_MINORDER
MAXORDER = FMMGEN_MAXORDER

cdef extern from "variant.hpp":
    cdef enum class FMMVariantKind "FMMVariantKind":
        Full
        Compressed
        Planar
    bool fmm_have_compressed() except +
    bool fmm_have_planar() except +
    void fmm_select(FMMVariantKind kind) except +

HAVE_COMPRESSED = fmm_have_compressed()
HAVE_PLANAR = fmm_have_planar()

# fmm_select's target is a single process-global pointer that
# compute_field_fmm/bh/exact read at CALL time, not at build_tree time (see
# the note in variant.hpp): an FMM instance's tree is built once with the
# coefficient array sizes for whichever variant was selected at that moment,
# but every later solve looks the operator functions up again through the
# same global. If a 2D (planar) and a 3D (full/compressed) instance are both
# alive in one process - exactly what DemagFMM auto-selecting by mesh
# dimensionality invites - whichever selected last would silently apply to
# both unless each instance re-selects its own variant immediately before
# every call. Both classes below do that rather than only once at
# construction.
cdef extern from "dim_shim.hpp":
    cdef cppclass Tree3:
        Tree3()
        void compute_field_fmm(double *F)
        void compute_field_exact(double *F)
        void compute_field_bh(double *F)
    cdef cppclass Tree2:
        Tree2()
        void compute_field_fmm(double *F)
        void compute_field_exact(double *F)
        void compute_field_bh(double *F)
    Tree3 build_tree_3d(double *pos, double *mu, size_t nparticles, size_t ncrit, size_t order, double theta)
    Tree2 build_tree_2d(double *pos, double *mu, size_t nparticles, size_t ncrit, size_t order, double theta)


cdef class FMM:
    """FMM/Barnes-Hut tree over general 3D (x, y, z) positions."""
    cdef public size_t n
    cdef public size_t ncrit
    cdef public double theta
    cdef public size_t order
    cdef public double [:, :] r
    cdef public double [:] mu
    cdef public double [:] mu_s
    cdef public double [:] Mu
    cdef public int calc_type
    cdef Tree3 tree
    cdef FMMVariantKind variant_kind

    def __cinit__(self, size_t n, size_t ncrit, double theta, size_t order, double [:, :] r, double [:] mu, double [:] mu_s, calc_type=0, compressed=True):
        if order < MINORDER or order >= MAXORDER:
            raise ValueError(
                f"Order needs to be in [{MINORDER}, {MAXORDER}), got {order}")
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

        self.variant_kind = FMMVariantKind.Compressed if compressed else FMMVariantKind.Full
        fmm_select(self.variant_kind)
        self.tree = build_tree_3d(&self.r[0, 0], &self.Mu[0], self.n, self.ncrit, self.order, self.theta)

    cdef _scale(self):
        for i in range(self.n):
            self.Mu[3*i + 0] = self.mu[3*i + 0] * self.mu_s[i]
            self.Mu[3*i + 1] = self.mu[3*i + 1] * self.mu_s[i]
            self.Mu[3*i + 2] = self.mu[3*i + 2] * self.mu_s[i]

    def compute_field(self, double [:] F):
        fmm_select(self.variant_kind)
        self._scale()
        if self.calc_type == 0:
            self.tree.compute_field_fmm(&F[0])
        elif self.calc_type == 1:
            self.tree.compute_field_bh(&F[0])

    def compute_field_exact(self, double [:] F):
        fmm_select(self.variant_kind)
        self._scale()
        self.tree.compute_field_exact(&F[0])


cdef class FMM2D:
    """FMM/Barnes-Hut tree over planar (x, y) positions, z implicitly 0.

    Source (dipole moment) vectors stay full 3D - a planar arrangement of
    dipoles can still point out of the plane - only positions are 2-wide.
    Always uses fmmgen's planar operator variant
    (generate_code(..., planar=True)), which exploits every particle having
    z=0; there is no full/compressed choice the way there is for FMM's
    general 3D case, since planar is the only sensible operator set once
    positions are genuinely 2D.
    """
    cdef public size_t n
    cdef public size_t ncrit
    cdef public double theta
    cdef public size_t order
    cdef public double [:, :] r
    cdef public double [:] mu
    cdef public double [:] mu_s
    cdef public double [:] Mu
    cdef public int calc_type
    cdef Tree2 tree

    def __cinit__(self, size_t n, size_t ncrit, double theta, size_t order, double [:, :] r, double [:] mu, double [:] mu_s, calc_type=0):
        if order < MINORDER or order >= MAXORDER:
            raise ValueError(
                f"Order needs to be in [{MINORDER}, {MAXORDER}), got {order}")
        self.calc_type = calc_type
        self.n = n
        self.ncrit = ncrit
        self.theta = theta
        self.order = order
        self.r = r
        self.mu = mu
        self.mu_s = mu_s
        self.Mu = np.zeros(3*self.n)

        fmm_select(FMMVariantKind.Planar)
        self.tree = build_tree_2d(&self.r[0, 0], &self.Mu[0], self.n, self.ncrit, self.order, self.theta)

    cdef _scale(self):
        for i in range(self.n):
            self.Mu[3*i + 0] = self.mu[3*i + 0] * self.mu_s[i]
            self.Mu[3*i + 1] = self.mu[3*i + 1] * self.mu_s[i]
            self.Mu[3*i + 2] = self.mu[3*i + 2] * self.mu_s[i]

    def compute_field(self, double [:] F):
        fmm_select(FMMVariantKind.Planar)
        self._scale()
        if self.calc_type == 0:
            self.tree.compute_field_fmm(&F[0])
        elif self.calc_type == 1:
            self.tree.compute_field_bh(&F[0])

    def compute_field_exact(self, double [:] F):
        fmm_select(FMMVariantKind.Planar)
        self._scale()
        self.tree.compute_field_exact(&F[0])
