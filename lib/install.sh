#!/bin/bash

SUNDIALS=sundials-2.5.0
FFTW=fftw-3.3.3

#============================================================================
#sundials
#============================================================================
if [ ! -e ${SUNDIALS}.tar.gz ]; then
  wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/${SUNDIALS}.tar.gz
fi;

if [ ! -e ${SUNDIALS} ]; then
  tar -xvzf ${SUNDIALS}.tar.gz
  cd ${SUNDIALS}
  ./configure --enable-shared
  make
  sudo make install
  cd ..
fi;

#============================================================================
# fftw3
#============================================================================

if [ ! -e ${FFTW}.tar.gz ]; then
  wget http://www.fftw.org/${FFTW}.tar.gz
fi;

if [ ! -e ${FFTW} ]; then
  tar -xvzf ${FFTW}.tar.gz
  cd ${FFTW}
  ./configure --enable-shared
  make
  sudo make install
  cd ..
fi;
