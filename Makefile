PROJECT_DIR = $(abspath .)
EXTENSIONS_DIR = ${PROJECT_DIR}/fidimag/extensions

#####################
# Cython Extensions #
#####################


build: extensions-directory
	python setup.py build_ext --build-lib ${EXTENSIONS_DIR}
	touch ${EXTENSIONS_DIR}/__init__.py
	rm -rf ${PROJECT_DIR}/build

extensions-directory:
	mkdir -p ${EXTENSIONS_DIR}

clean:
	rm -rf ${EXTENSIONS_DIR}

#########
# Tests #
#########

create-dirs:
	mkdir -p test-reports/junit

test: create-dirs
	py.test -v --junitxml=$(PROJECT_DIR)/test-reports/junit/test-pytest.xml

test-without-run-oommf: create-dirs
	py.test -v -s --junitxml=$(PROJECT_DIR)/test-reports/junit/test-pytest.xml -m "not run_oommf"

test-basic:
	cd tests && py.test -v

test-ipynb: create-dirs
	cd doc/ipynb && py.test . -v --ipynb --sanitize-with sanitize_file --junitxml=$(PROJECT_DIR)/test-reports/junit/test-ipynb-pytest.xml

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
