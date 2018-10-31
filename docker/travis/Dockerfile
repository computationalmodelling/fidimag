FROM ubuntu:16.04

# To build this image `docker build -t fidimag .`
# Then you can drop into a live bash shell with `docker run -it fidimag`.
ENTRYPOINT ["/bin/bash"]  
SHELL ["/bin/bash", "-c"] 

RUN apt-get update -y
RUN apt-get install -y build-essential cmake cython3 python3-dev python3-pip \
    python3-psutil python3-pytest python3-pytest-cov \
    python3-matplotlib python3-numpy python3-scipy \
    liblapack-dev libopenblas-dev \
    wget
# next line for codecov target
RUN apt-get install -y curl git

# ----------
# hack to fix sudden breakage in CI
# (fix from https://github.com/getsentry/sentry/issues/3143)
# (first occurance of fail at https://travis-ci.org/computationalmodelling/fidimag/builds/319708056?utm_source=github_status&utm_medium=notification)
# (which is part of this pull request: https://github.com/computationalmodelling/fidimag/pull/106)
RUN pip3 install --upgrade setuptools==20.4
# ----------

# ----------
# hack to fix sudden breakage in CI
# (fix from https://github.com/getsentry/sentry/issues/3143)
# (first occurance of fail at https://travis-ci.org/computationalmodelling/fidimag/builds/319708056?utm_source=github_status&utm_medium=notification)
# (which is part of this pull request: https://github.com/computationalmodelling/fidimag/pull/106)
RUN pip3 install --upgrade setuptools==20.4
# ----------

RUN pip3 install ipywidgets nbval pyvtk six pytest -U

WORKDIR /fidimag
ADD . /fidimag
RUN ./bin/install-sundials.sh
RUN ./bin/install-fftw.sh
RUN make build

ENV PYTHONPATH=/fidimag \
    LD_LIBRARY_PATH=/fidimag/local/lib LD_RUN_PATH=/fidimag/local/lib \
    OMP_NUM_THREADS=1 MPLBACKEND=Agg QT_API=pyqt
