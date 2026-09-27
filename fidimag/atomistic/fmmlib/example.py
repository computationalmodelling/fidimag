import fmmgen

# Always run this script from inside fmmlib/ (not a temp/build directory
# copied back afterwards): the generated `#include` uses whatever path is
# in front of "operators" below, so running it elsewhere leaves an absolute
# path in line 1, and CMake compiles the checked-in file in place, expecting
# `#include "operators.h"`.
#
# The order-dispatch wrappers at the end of the file used to paste the
# declaration's `__restrict` into the *call*, e.g.
# `S2M_2(x, y, z, __restrict S, __restrict M);`, which is not valid C++.
# This is fixed in fmmgen's writer (rebuilds each dispatch call from the
# declaration's parameter names, and uses an `FMMGEN_RESTRICT` macro so the
# same generated code is valid whether compiled as C or C++) -- but that fix
# is, as of this regeneration, an uncommitted patch on the local `fmmgen`
# checkout's `fidimag-2d-vendor` branch, not yet upstream. Regenerating
# against a fmmgen checkout without that patch will reintroduce the bug; if
# so, the call-site fixup below is the manual workaround used previously:
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
