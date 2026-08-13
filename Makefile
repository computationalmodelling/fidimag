PROJECT_DIR = $(abspath .)
EXTENSIONS_DIR = ${PROJECT_DIR}/fidimag/extensions
USER_EXTENSIONS_DIR = ${PROJECT_DIR}/fidimag/extensions/user
# Use uv run to execute in the managed virtual environment
PYTHON = uv run python
PYTEST = uv run pytest

.DEFAULT_GOAL := help

help:
	@echo "Fidimag Makefile - Common Tasks"
	@echo "================================"
	@echo ""
	@echo "Build:"
	@echo "  make build          Build fidimag using uv (recommended)"
	@echo "  make build-old      Build using old setup.py (deprecated)"
	@echo "  make clean          Clean all build artifacts"
	@echo ""
	@echo "Testing:"
	@echo "  make test           Run quick tests (not slow, not oommf)"
	@echo "  make test-all       Run all tests"
	@echo "  make test-basic     Run basic tests"
	@echo "  make test-oommf     Run OOMMF tests only"
	@echo "  make test-ipynb     Run notebook tests"
	@echo "  make test-clean     Clean test output files"
	@echo ""
	@echo "Documentation:"
	@echo "  make doc            Build all documentation"
	@echo "  make doc-html       Build HTML documentation"
	@echo "  make doc-clean      Clean documentation build"
	@echo ""
	@echo "Docker:"
	@echo "  make docker         Build Docker image"
	@echo "  make test-docker    Run tests in Docker"
	@echo ""
	@echo "Note: All Python/pytest commands use 'uv run' to execute in the managed environment"
	@echo ""

#####################
# Cython Extensions #
#####################

build:
	uv sync --reinstall-package fidimag

build-old:
	# Old setup.py-based build (deprecated)
	${PYTHON} setup.py build_ext --inplace -j2

clean:
	bash clean_cython_files.sh

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
