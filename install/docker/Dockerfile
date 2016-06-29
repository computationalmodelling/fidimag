FROM ubuntu:14.04

# packages we need to run fidimag
RUN apt-get -y update
RUN apt-get -y install python-numpy python-dev python-scipy
RUN apt-get -y install python-pytest python-pyvtk ipython python-matplotlib mayavi2
# standard tools for compilation
RUN apt-get -y install wget make git

# where to install source
ENV FIDIMAG_HOME /home/fidimag

RUN mkdir -p $FIDIMAG_HOME
WORKDIR $FIDIMAG_HOME
RUN git clone https://github.com/computationalmodelling/fidimag.git
WORKDIR $FIDIMAG_HOME/fidimag/bin

# install third party libraries from source
RUN bash install-ubuntu-packages.sh
RUN bash install-fftw.sh
RUN bash install-sundials-2.5.sh

# for pip
RUN apt-get -y install python-pip
# install cython
RUN pip install cython --upgrade
WORKDIR $FIDIMAG_HOME/fidimag

# compile fidimag
RUN make
env PYTHONPATH=$FIDIMAG_HOME/fidimag
env LD_LIBRARY_PATH=$FIDIMAG_HOME/fidimag/local/lib
WORKDIR $FIDIMAG_HOME/fidimag/tests

# check that tests run okay
RUN py.test -v


# install Jupyter
RUN pip install jupyter


# Set up user so that we do not run as root
RUN useradd -m -s /bin/bash -G sudo fidimag && \
    echo "fidimag:docker" | chpasswd && \
    echo "fidimag ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
RUN chown -R fidimag $FIDIMAG_HOME

USER fidimag
RUN touch $FIDIMAG_HOME/.sudo_as_admin_successful

# Print something nice on entry.
#COPY WELCOME $FIDIMAG_HOME/WELCOME


WORKDIR $FIDIMAG_HOME
CMD ["/bin/bash","-i"]
