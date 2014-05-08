all:
	python setup.py build_ext --inplace


iridis:
	python setup_iridis.py build_ext --inplace

	#rebuild clib.so and cvode.so using icc
	icc -pthread -shared build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/clib.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/anis.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/demag.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/dmi.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/exch.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/llg.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/sllg.o build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/stt.o -L/home/ww1g11/Softwares/pccp/libs/lib -L/local/software/intel/2013.4.183/composer_xe_2013.4.183/compiler/lib/intel64 -Wl,-R/home/ww1g11/Softwares/pccp/libs/lib -Wl,-R/local/software/intel/2013.4.183/composer_xe_2013.4.183/compiler/lib/intel64 -lm -lfftw3 -lsundials_cvodes -lsundials_nvecserial -liomp5 -o /home/ww1g11/Softwares/pccp/clib.so

	icc -pthread -shared build/temp.linux-x86_64-2.7/home/ww1g11/Softwares/pccp/cp/sundials/cvode.o -L/home/ww1g11/Softwares/pccp/libs/lib -L/local/software/intel/2013.4.183/composer_xe_2013.4.183/compiler/lib/intel64 -Wl,-R/home/ww1g11/Softwares/pccp/libs/lib -Wl,-R/local/software/intel/2013.4.183/composer_xe_2013.4.183/compiler/lib/intel64 -lm -lfftw3 -lsundials_cvodes -lsundials_nvecserial -liomp5 -o /home/ww1g11/Softwares/pccp/cvode.so
	
