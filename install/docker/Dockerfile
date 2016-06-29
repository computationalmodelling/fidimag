FROM ubuntu:14.04
# Do these first as I expect them to go wrong first ;)
# required to compile fidimag
RUN apt-get -y update
RUN apt-get -y install python-numpy python-dev python-scipy
RUN apt-get -y install python-pytest python-pyvtk ipython python-matplotlib mayavi2
# standard tools for compilation
RUN apt-get -y install wget make git

ENV FIDIMAG_HOME /home/fidimag

RUN mkdir -p $FIDIMAG_HOME
WORKDIR $FIDIMAG_HOME
RUN git clone https://github.com/computationalmodelling/fidimag.git
RUN pwd
RUN cd fidimag/bin && pwd
WORKDIR $FIDIMAG_HOME/fidimag/bin
RUN pwd

RUN bash install-ubuntu-packages.sh
RUN bash install-fftw.sh
RUN bash install-sundials-2.5.sh
# for pip
RUN apt-get -y install python-pip
RUN pip install cython --upgrade
WORKDIR $FIDIMAG_HOME/fidimag
RUN make
env PYTHONPATH=$FIDIMAG_HOME/fidimag
env LD_LIBRARY_PATH=$FIDIMAG_HOME/fidimag/local/lib
WORKDIR $FIDIMAG_HOME/fidimag/tests
# check that tests run okay
RUN py.test -v

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
#USER root
#ENTRYPOINT ["/sbin/my_init","--quiet","--","/sbin/setuser","fidimag","/bin/bash","-l","-c"]
#ENTRYPOINT ["/sbin/setuser","fidimag","/bin/bash","-l","-c"]
CMD ["/bin/bash","-i"]

