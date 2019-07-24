import fmmgen

order = 15

fmmgen.generate_code(order, "operators", CSE=True, cython_wrapper=False, potential=False, field=True, minpow=5)
