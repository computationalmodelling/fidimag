#!/bin/bash

SUNDIALS=sundials-2.5.0

#============================================================================
#sundials
#============================================================================
if [ ! -e ${SUNDIALS}.tar.gz ]; then
  wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/sundials-2.5.0.tar.gz
  tar -xvzf ${SUNDIALS}.tar.gz
  cd ${SUNDIALS}
  ./configure
  make
  sudo make install
fi;



#the same as fftw. will add later ....