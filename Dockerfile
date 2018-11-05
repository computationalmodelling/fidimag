FROM ubuntu:16.04

RUN apt -y update 
RUN apt install -y git python3 python3-pip gcc psutils cmake wget make
RUN apt install -y gfortran libblas-dev liblapack-dev python3-tk sudo fonts-lato
RUN pip3 install cython matplotlib pytest scipy psutil pyvtk ipywidgets -U
RUN pip3 install --no-cache-dir notebook

RUN ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /usr/local
RUN git clone https://github.com/computationalmodelling/fidimag.git
WORKDIR /usr/local/fidimag
# Work with stable release
RUN git checkout tags/v2.9
# Install CVODE and FFTW libraries
WORKDIR /usr/local/fidimag/bin
RUN bash install-fftw.sh
RUN bash install-sundials.sh

ENV PYTHONPATH="/usr/local/fidimag:$PYTHONPATH"
ENV LD_LIBRARY_PATH="/usr/local/fidimag/local/lib:$LD_LIBRARY_PATH"

WORKDIR /usr/local/fidimag
RUN python3 setup.py build_ext --inplace
RUN python3 -c "import matplotlib"
# Headless Matplotlib:
ENV MPLBACKEND=Agg

# Headless Matplotlib:
ENV MPLBACKEND=Agg

# Set threads for OpenMP:
ENV OMP_NUM_THREADS=2
# WORKDIR /io

# User to make Binder happy
ENV NB_USER magnetism
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

WORKDIR /home/${USER}/magnetism/doc/ipynb
