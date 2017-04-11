FROM ubuntu:14.04

# packages we need to run fidimag
RUN apt-get -y update
RUN apt-get -y install python3-numpy python3-dev python3-scipy cmake
RUN apt-get -y install python3-pytest ipython3 python3-matplotlib python3-pip
# standard tools for compilation
RUN apt-get -y install wget make git

# where to install source
ENV FIDIMAG_HOME /home/fidimag

RUN mkdir -p $FIDIMAG_HOME
WORKDIR $FIDIMAG_HOME
RUN git clone https://github.com/computationalmodelling/fidimag.git
WORKDIR $FIDIMAG_HOME/fidimag/bin

# install third party libraries from source
RUN bash install-fftw.sh
RUN bash install-sundials.sh

# for pip
RUN python3 -m pip install --user --upgrade setuptools pip
# install pyvtk
RUN python3 -m pip install pyvtk
# install cython
RUN python3 -m pip install cython --upgrade
WORKDIR $FIDIMAG_HOME/fidimag
RUN apt-get -y install libatlas-base-dev libatlas3-base libatlas3gf-base
# compile fidimag
RUN python3 setup.py build_ext --inplace
env PYTHONPATH=$FIDIMAG_HOME/fidimag
env LD_LIBRARY_PATH=$FIDIMAG_HOME/fidimag/local/lib
WORKDIR $FIDIMAG_HOME/fidimag/tests

# check that tests run okay
RUN py.test-3 -v


# install Jupyter, port exposing doesn't work yet
#RUN pip install jupyter

# expose jupyter port - not working yet
#EXPOSE 8888 8888


# Set up user so that we do not run as root
RUN useradd -m -s /bin/bash -G sudo fidimag && \
    echo "fidimag:docker" | chpasswd && \
    echo "fidimag ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
RUN chown -R fidimag $FIDIMAG_HOME

# For bind mounts from host
RUN mkdir /io
RUN chown -R fidimag /io

USER fidimag
RUN touch $FIDIMAG_HOME/.sudo_as_admin_successful

RUN cd fidimag
RUN py.test -v
# Print something nice on entry.
#COPY WELCOME $FIDIMAG_HOME/WELCOME

WORKDIR /io
CMD ["/bin/bash","-i"]
