.PHONY: build

# build cython files

PROJECT_DIR = $(abspath .)

build:
	python setup.py build_ext --inplace

create-dirs:
	mkdir -p test-reports/junit

clean:
	rm *.so

test: create-dirs
	py.test -v --junitxml=$(PROJECT_DIR)/test-reports/junit/test-pytest.xml

test-micro:
	cd micro/tests && py.test -v

test-basic:
	cd tests && py.test -v

test-ipynb: create-dirs
	cd doc/ipynb && py.test -v --ipynb . --junitxml=$(PROJECT_DIR)/test-reports/junit/test-ipynb-pytest.xml

# Building documentation
doc: doc-html doc-latexpdf doc-singlehtml

doc-clean:
	make -C doc clean

doc-%:
	@echo $*
	make -C doc $*
