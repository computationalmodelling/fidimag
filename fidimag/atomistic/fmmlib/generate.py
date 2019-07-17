import fmmgen

order = 13

fmmgen.generate_code(order, "operators", CSE=True, cython_wrapper=False, potential=False, field=True)

