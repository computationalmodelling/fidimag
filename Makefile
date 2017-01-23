PROJECT_DIR = $(abspath .)
EXTENSIONS_DIR = ${PROJECT_DIR}/fidimag/extensions
PYTHON = python

#####################
# Cython Extensions #
#####################


build:
	${PYTHON} setup.py build_ext --inplace

clean:
	rm -rf ${EXTENSIONS_DIR}/*
	touch ${EXTENSIONS_DIR}/__init__.py

#########
# Tests #
#########

# Quick tests, also not using OOMMF tests
test:
	cd tests && py.test -v -m "not slow and not run_oommf"

test-clean:
	rm -rf neb*.ndt
	rm -rf relax_npys
	rm -rf relax.txt
	rm -rf skx_number_*
	rm -rf skx_num_*
	rm -rf vtks
	rm -rf spin.txt
	rm -rf unnamed_npys
	rm -rf unnamed.txt
	rm -rf 2dpbc_tensors.npz
	rm -rf 10spin.txt
	rm -rf dyn_npys
	rm -rf dyn_spin.txt
	rm -rf dyn.txt
	rm -rf dyn_vtks
	rm -rf m0.npy
	rm -rf skx_num.txt
	rm -rf test_energy.txt
	rm -rf unnamed_vtks
	rm -rf dw_cobalt.txt
	rm -rf mu_s.npy
	rm -rf npys
	rm -rf relax_vtks
	rm -rf m1.npy

test2:
	# like test, but run also outside the 'tests' directory.
	# Doesn't work on Hans laptop.
	py.test -v -m "not slow and not run_oommf"

test-all: create-dirs
	py.test -v --junitxml=$(PROJECT_DIR)/test-reports/junit/test-pytest.xml

test-without-run-oommf: create-dirs
	py.test -v -m "not run_oommf" --cov=fidimag --cov-report=html --junitxml=$(PROJECT_DIR)/test-reports/junit/test-pytest.xml

test-basic:
	cd tests && py.test -v

# Convenience name for commonly used quick running of tests
tq:
	$(error This target 'tq' has been removed, please update the code calling this)

test-quick:
	$(error This target 'test-quick' has been removed, please update the code calling this)


test-ipynb: create-dirs
	cd doc/ipynb && py.test . -v --current-env --nbval --sanitize-with sanitize_file --junitxml=$(PROJECT_DIR)/test-reports/junit/test-ipynb-pytest.xml

test-oommf:
	py.test -v -m "oommf"

create-dirs:
	mkdir -p test-reports/junit


#################
# Documentation #
#################

doc: doc-html doc-latexpdf doc-singlehtml

doc-clean:
	make -C doc clean

doc-%:
	@echo $*
	make -C doc $*

.PHONY: extensions-directory build clean create-dirs test test-basic test-ipynb doc doc-clean doc-html doc-latexpdf doc-singlehtml
