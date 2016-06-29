FROM ubuntu:14.04
# Do these first as I expect them to go wrong first ;)
# required to compile fidimag
RUN apt-get -y update
RUN apt-get -y install python-numpy python-dev python-scipy
RUN apt-get -y install python-pytest python-pyvtk ipython python-matplotlib mayavi2
# standard tools for compilation
RUN apt-get -y install wget make git

RUN mkdir -p /mnt/install && cd /mnt/install
RUN git clone https://github.com/computationalmodelling/fidimag.git
RUN cd fidimag/bin

RUN bash install-ubuntu-packages.sh
RUN bash install-fftw.sh
RUN bash install install-sundials-2.5.sh
# for pip
RUN apt-get -y install python-pip
RUN pip install cython --upgrade
RUN cd ..
RUN make
env PYTHONPATH=/mnt/install/fidimag
env LD_LIBRARY_PATH=/mnt/install/fidimag/local/lib
RUN cd tests
CMD cd /mnt/install/fidimag/tests && py.test -v


