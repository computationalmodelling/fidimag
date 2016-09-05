# required to compile fidimag
deps_compilation="python-pip python-numpy python-dev python-scipy cmake"

# required for tests and running fidimag
deps_live="python-pytest python-pyvtk ipython python-matplotlib"

apt-get install $deps_compilation $deps_live
