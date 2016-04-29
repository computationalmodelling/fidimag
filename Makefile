PROJECT_DIR = $(abspath .)
EXTENSIONS_DIR = ${PROJECT_DIR}/fidimag/extensions

#####################
# Cython Extensions #
#####################


build:
	python setup.py build_ext --inplace

clean:
	rm -rf ${EXTENSIONS_DIR}/*
	touch ${EXTENSIONS_DIR}/__init__.py

#########
# Tests #
#########

# Quick tests, also not using OOMMF tests
test:
	cd tests && py.test -v -m "not slow and not run_oommf"

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
	cd doc/ipynb && py.test . -v --ipynb --sanitize-with sanitize_file --junitxml=$(PROJECT_DIR)/test-reports/junit/test-ipynb-pytest.xml

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
