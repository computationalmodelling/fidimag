import fmmgen

order = 11

fmmgen.generate_code(order, "operators", CSE=True, cython_wrapper=False, potential=False, field=True, minpower=5)
