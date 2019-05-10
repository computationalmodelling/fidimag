PROJECT_DIR = $(abspath .)
EXTENSIONS_DIR = ${PROJECT_DIR}/fidimag/extensions
PYTHON = python3
PYTEST = ${PYTHON} -m pytest

CC=gcc
CXX=g++

FFTW_INC=local/include
FFTW_LIB=local/lib

SUNDIALS_INC=local/include
SUNDIALS_LIB=local/lib

CPPFLAGS = -fPIC -g -Iinclude/ \
           -I${FFTW_INC} -L${FFTW_LIB} \
           -I${SUNDIALS_INC} -L${SUNDIALS_LIB} \
           -Inative/include \
           -lm \
           -lfftw3_omp -lfftw3 \
           -lsundials_cvodes -lsundials_nvecserial -lsundials_nvecopenmp \
           -lblas -llapack \
           -fopenmp

LDFLAGS = -shared
SOURCES = $(shell echo native/src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
LIBRARY = local/lib/libfidimag.so

########
# Build
########


all: $(LIBRARY) build


$(LIBRARY) : $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)


build: $(LIBRARY)
	CC=${CC} CXX=${CXX} CPPFLAGS="${CPPFLAGS}" ${PYTHON} setup.py build_ext --inplace


clean:
	rm -rf ${EXTENSIONS_DIR}/*
	touch ${EXTENSIONS_DIR}/__init__.py
	rm -rf build
	rm -rf $(OBJECTS) $(TARGET) *.dSYM
	find fidimag/ "*.cpp" -exec echo {} \;


docker:
	docker build -t fidimag -f ./docker/travis/Dockerfile .
	docker run -ti -d --name fidimag fidimag

#########
# Tests #
#########

test-docker:
	docker exec fidimag make test-basic
	#docker exec fidimag make test-without-run-oommf
	#docker exec fidimag make test-ipynb

ipynb-docker:
	while sleep 9m; do echo "===[ Still Running ]===\n"; done &
	docker exec fidimag make test-ipynb
	kill %1

travis-test: test-docker
	docker exec fidimag make codecov
	docker exec fidimag make test-ipynb

codecov:
	bash <(curl -s https://codecov.io/bash)
codecov: SHELL:= /bin/bash
# or the recipe fails with /bin/sh complaining
# `/bin/sh: 1: Syntax error: "(" unexpected`.

# Quick tests, also not using OOMMF tests
test:
	cd tests && ${PYTEST} -v -m "not slow and not run_oommf"

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
	${PYTEST} -v -m "not slow and not run_oommf"

test-all: create-dirs
	${PYTEST} -v --junitxml=$(PROJECT_DIR)/test-reports/junit/test-pytest.xml

test-without-run-oommf: create-dirs
	${PYTEST} -v -m "not run_oommf" --cov=fidimag --cov-report=html --junitxml=$(PROJECT_DIR)/test-reports/junit/test-pytest.xml

test-basic:
	cd tests && ${PYTEST} -v

# Convenience name for commonly used quick running of tests
tq:
	$(error This target 'tq' has been removed, please update the code calling this)

test-quick:
	$(error This target 'test-quick' has been removed, please update the code calling this)


test-ipynb: create-dirs
	cd doc/ipynb && ${PYTEST} . -v --current-env --nbval --sanitize-with sanitize_file --junitxml=$(PROJECT_DIR)/test-reports/junit/test-ipynb-pytest.xml

test-oommf:
	${PYTEST} -v -m "oommf"

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

.PHONY: extensions-directory build clean create-dirs test test-basic test-ipynb doc doc-clean doc-html doc-latexpdf doc-singlehtml docker test-docker travis codecov
