FROM jupyter/scipy-notebook
# where to install source
USER root
RUN mkdir -p /io
RUN chown $NB_USER /io
RUN apt update && apt install -y cmake
USER $NB_USER
ENV FIDIMAG_DIR /io
WORKDIR /io
RUN git clone https://github.com/rpep/fidimag.git
WORKDIR /io/fidimag
# install third party libraries from source
RUN bash bin/install-fftw.sh
RUN bash bin/install-sundials.sh

# install pyvtk
RUN pip install pyvtk
# install cython
RUN pip install cython --upgrade

# compile fidimag
RUN python3 setup.py build_ext --inplace
RUN pip install psutil
ENV PYTHONPATH=$FIDIMAG_DIR
ENV LD_LIBRARY_PATH=$FIDIMAG_DIR/local/lib
WORKDIR $FIDIMAG_DIR/tests

# https://github.com/conda-forge/matplotlib-feedstock/issues/36
RUN conda install --quiet --yes icu

# check that tests run okay
RUN conda install --quiet --yes pytest

# /io will be mounted from the host system
WORKDIR /io

