.PHONY: build

# build cython files

PROJECT_DIR = $(abspath .)

build:
	exit 17 # Testing a build error or jenkins. Try another commit if you encounter this.
	python setup.py build_ext --inplace

create-dirs:
	mkdir -p test-reports/junit

clean:
	rm *.so

test: create-dirs
	py.test -v --junitxml=$(PROJECT_DIR)/test-reports/junit/TEST_pytest.xml

test-micro:
	cd micro/tests && py.test -v

test-basic:
	cd tests && py.test -v

# Building documentation
doc: doc-html doc-latexpdf doc-singlehtml

doc-clean:
	make -C doc clean

doc-%:
	@echo $*
	make -C doc $*
