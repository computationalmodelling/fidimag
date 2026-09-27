#pragma once
#include <cstddef>

/*! \brief Runtime selection between the uncompressed, harmonic-compressed and
    planar (2D-plane) operator sets.

    fmmgen can emit up to three operator families into operators.h: the
    unsuffixed ones (Nterms(order) coefficients per expansion), c-suffixed
    ones when generated with compress=True (the trace-free (order+1)^2), and
    xy-suffixed ones when generated with planar=True (sources, and for the
    local array's retained set, targets confined to z=0; moment vectors stay
    full 3D). The three compute the same field under their respective
    validity conditions -- compressed always, planar only when the tree's
    particles genuinely have z=0 -- so which one runs is purely a performance
    question, worth flipping without rebuilding.

    One binary holding all three, selected through this indirection, means an
    A/B/C comparison changes only the operators. The indirection costs one
    indirect call per operator invocation; the target is invariant for the
    whole run so it predicts perfectly, and each call does hundreds of flops.

    Selection must happen before build_tree, which sizes the coefficient
    arrays from fmm->msize()/lsize().
*/
struct FMMVariant {
  void (*s2m)(double, double, double, double *, double *, int);
  void (*m2m)(double, double, double, double *, double *, int);
  void (*m2l)(double, double, double, double *, double *, int);
  void (*l2l)(double, double, double, double *, double *, int);
  void (*l2p)(double, double, double, double *, double *, int);
  void (*m2p)(double, double, double, double *, double *, int);
  size_t (*msize)(size_t order);   //!< coefficients per multipole expansion
  size_t (*lsize)(size_t order);   //!< coefficients per local expansion
  const char *name;
};

enum class FMMVariantKind { Full, Compressed, Planar };

//! The selected operator set. Defaults to uncompressed.
extern const FMMVariant *fmm;

//! True if operators.h was generated with compress=True.
bool fmm_have_compressed();

//! True if operators.h was generated with planar=True.
bool fmm_have_planar();

/*! Select the operator set. Throws if the requested set was not generated --
    silently falling back would make a benchmark quietly compare a thing
    against itself. Planar operators are only physically valid if every
    particle's z coordinate is actually 0 (see tree.hpp's D=2 instantiation);
    this function has no way to check that and does not try to. */
void fmm_select(FMMVariantKind kind);
