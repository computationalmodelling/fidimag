import fmmgen

order = 8
source_order = 1
cse = True
atomic = True
precision='double'
fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     cython=False,
                     harmonic_derivs=True,
                     potential=False,
                     field=True,
                     source_order=source_order,
                     atomic=atomic)
