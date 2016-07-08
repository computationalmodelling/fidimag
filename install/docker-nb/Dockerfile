FROM jupyter/scipy-notebook

# where to install source
ENV FIDIMAG_DIR $HOME/work/fidimag

RUN git clone https://github.com/computationalmodelling/fidimag.git
WORKDIR $FIDIMAG_DIR

# install third party libraries from source
RUN bash bin/install-fftw.sh
RUN bash bin/install-sundials-2.5.sh

# install pyvtk
RUN pip install pyvtk
# install cython
RUN pip install cython --upgrade

# compile fidimag
RUN python3 setup.py build_ext --inplace
ENV PYTHONPATH=$FIDIMAG_DIR
ENV LD_LIBRARY_PATH=$FIDIMAG_DIR/local/lib
WORKDIR $FIDIMAG_DIR/tests

# https://github.com/conda-forge/matplotlib-feedstock/issues/36
RUN conda install --quiet --yes icu

# check that tests run okay
RUN conda install --quiet --yes pytest
RUN py.test -v


WORKDIR $HOME/work
