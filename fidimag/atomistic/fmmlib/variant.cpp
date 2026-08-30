#include "variant.hpp"
#include "operators.h"
#include "utils.hpp"
#include <stdexcept>

static size_t full_msize(size_t order) { return Msize(order, FMMGEN_SOURCEORDER); }
static size_t full_lsize(size_t order) { return Lsize(order, FMMGEN_SOURCEORDER); }

static const FMMVariant kFull = {
    S2M, M2M, M2L, L2L, L2P, M2P, full_msize, full_lsize, "uncompressed",
};

#ifdef FMMGEN_COMPRESSED
// The generated tables, rather than a local (order+1)*(order+1): the header is
// the single source of truth for sizes the caller cannot otherwise derive.
static size_t comp_msize(size_t order) { return FMMGEN_MULTIPOLESIZE[order]; }
static size_t comp_lsize(size_t order) { return FMMGEN_LOCALSIZE[order]; }

static const FMMVariant kCompressed = {
    S2Mc, M2Mc, M2Lc, L2Lc, L2Pc, M2Pc, comp_msize, comp_lsize, "harmonic-compressed",
};
#endif

const FMMVariant *fmm = &kFull;

bool fmm_have_compressed() {
#ifdef FMMGEN_COMPRESSED
  return true;
#else
  return false;
#endif
}

void fmm_select(bool compressed) {
#ifdef FMMGEN_COMPRESSED
  fmm = compressed ? &kCompressed : &kFull;
#else
  if (compressed) {
    throw std::runtime_error(
        "compressed operators were not generated: set compress=True in "
        "example/example.py and regenerate");
  }
  fmm = &kFull;
#endif
}
