#pragma once
#include <cstddef>

/*! \brief Runtime selection between the uncompressed and harmonic-compressed
    operator sets.

    fmmgen emits both into operators.h when generated with compress=True: the
    unsuffixed operators, which take Nterms(order) coefficients per expansion,
    and c-suffixed ones taking the trace-free (order+1)^2. The arithmetic they
    perform is different but the field they produce is the same, so which one
    the solve uses is purely a performance question -- and therefore one worth
    being able to flip without rebuilding.

    Two binaries and a compile-time switch would have been simpler, but the
    expected end-to-end effect is on the order of 1.2x, and code layout differs
    enough between separately-linked binaries to muddy a difference that small.
    Holding both in one image means an A/B comparison changes the operators and
    nothing else.

    The indirection costs one indirect call per operator invocation. The target
    is invariant for the whole run so it predicts perfectly, and each call does
    hundreds of flops; measure it rather than assume it (run --compress=0
    against a build without this layer).

    Selection must happen before build_tree, which sizes the coefficient arrays.
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

//! The selected operator set. Defaults to uncompressed.
extern const FMMVariant *fmm;

//! True if operators.h was generated with compress=True.
bool fmm_have_compressed();

/*! Select the operator set. Throws if the compressed set was requested but
    the operators were not generated with compress=True -- silently falling
    back would make a benchmark quietly compare a thing against itself. */
void fmm_select(bool compressed);
