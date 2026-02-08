FROM ubuntu:22.04

# To build this image `docker build -t fidimag .`
# Then you can drop into a live bash shell with `docker run -it fidimag`.
ENTRYPOINT ["/bin/bash"]
SHELL ["/bin/bash", "-c"]

ARG DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update -y -qq
RUN apt-get install -y -qq build-essential cmake gcc g++ gfortran \
    liblapack-dev libopenblas-dev \
    wget curl git ca-certificates

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    ln -s /root/.cargo/bin/uv /usr/local/bin/uv
ENV PATH="/root/.cargo/bin:${PATH}"

# Copy fidimag source
WORKDIR /fidimag
ADD . /fidimag

# Install SUNDIALS and FFTW to ./local/
RUN ./bin/install-sundials.sh
RUN ./bin/install-fftw.sh

# Build fidimag using uv with dev dependencies
RUN uv sync --frozen --all-extras

ENV LD_LIBRARY_PATH=/fidimag/local/lib \
    LD_RUN_PATH=/fidimag/local/lib \
    OMP_NUM_THREADS=1 \
    MPLBACKEND=Agg \
    QT_API=pyqt

# Run tests during build to verify
RUN uv run pytest tests/ -v -m "not slow and not run_oommf"
