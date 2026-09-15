import fmmgen

# Two fixups are needed on the operators.cpp this script writes, until both
# are fixed upstream in fmmgen's writer:
#
#   1. The generated `#include` uses whatever path was passed as the output
#      prefix, so running this outside fmmlib/ leaves an absolute path in
#      line 1. It must read `#include "operators.h"` for the CMake build,
#      which compiles the checked-in file in place.
#
#   2. The order-dispatch wrappers at the end of the file paste the
#      declaration's `__restrict` into the *call*, e.g.
#      `S2M_2(x, y, z, __restrict S, __restrict M);`, which is not valid
#      C++. Strip `__restrict` from the call sites only (the parameter
#      declarations it precedes a `*` in are correct and must stay):
#
#          perl -i -pe 's/(?<!\* )__restrict //g' operators.cpp
#

# generate_code's order argument is an EXCLUSIVE upper bound (FMMGEN_MAXORDER
# in the generated operators.h): passing order=13 generates orders 2..12
# inclusive, not order=12. The previous order=8 here (generating only 2..7)
# was this same off-by-one, compounded by a separate one in fmm.pyx's order
# validation that let requests for the ungenerated order=8 through silently
# instead of raising - see DemagFMM's order assertion in demag.py.
order = 13
source_order = 1
cse = True
atomic = True
precision = 'double'
fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     cython=False,
                     harmonic_derivs=True,
                     compress=True,
                     planar=True,
                     potential=False,
                     field=True,
                     source_order=source_order,
                     atomic=atomic,
                     # Horner-form the coordinate-polynomial arrays (S2M,
                     # M2M, L2L, L2P, M2P, and M2L's D) before CSE. Measured
                     # 1.6-3.7x faster at runtime on the compiled operators
                     # (M2M/L2L/L2P/M2P; M2L's own contraction is
                     # deliberately excluded, see fmmgen.writer's `horner`
                     # docstring) for ~20 minutes of extra one-time
                     # generation time at this order -- worth it since this
                     # script runs once to produce the checked-in
                     # operators.h/operators.cpp, not on every build.
                     horner=True,
                     language='c++')
