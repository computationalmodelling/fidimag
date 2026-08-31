import fmmgen

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
                     potential=False,
                     field=True,
                     source_order=source_order,
                     atomic=atomic,
                     language='c++')
