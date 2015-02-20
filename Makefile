# build cython files

PROJECT_DIR = $(abspath .)

build:
	exit 17 # Testing a build error or jenkins. Try another commit if you encounter this.
	python setup.py build_ext --inplace

create-dirs:
	mkdir -p test-reports/junit

iridis:
	# we will write a new one later when we have time
	python setup_iridis.py build_ext --inplace

	# rebuild clib.so and cvode.so using icc
	icc -pthread -shared build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/clib.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/anis.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/demag.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/dmi.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/exch.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/llg.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/sllg.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/stt.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/util.o -L/home/ww1g11/Softwares/pccp/libs/lib -L/local/software/intel/2013.4.183/composer_xe_2013.4.183/compiler/lib/intel64 -Wl,-R/home/ww1g11/Softwares/pccp/libs/lib -Wl,-R/local/software/intel/2013.4.183/composer_xe_2013.4.183/compiler/lib/intel64 -lm -lfftw3 -lsundials_cvodes -lsundials_nvecserial -liomp5 -o /home/ww1g11/Softwares/pccp/clib.so

	icc -pthread -shared build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/sundials/cvode.o -L/home/ww1g11/Softwares/pccp/libs/lib -L/local/software/intel/2013.4.183/composer_xe_2013.4.183/compiler/lib/intel64 -Wl,-R/home/ww1g11/Softwares/pccp/libs/lib -Wl,-R/local/software/intel/2013.4.183/composer_xe_2013.4.183/compiler/lib/intel64 -lm -lfftw3 -lsundials_cvodes -lsundials_nvecserial -liomp5 -o /home/ww1g11/Softwares/pccp/cvode.so

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
