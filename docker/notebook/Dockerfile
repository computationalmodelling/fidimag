FROM jupyter/scipy-notebook

USER root

# Install system dependencies
RUN apt update && apt install -y cmake gcc g++ build-essential curl ca-certificates

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    ln -s /root/.cargo/bin/uv /usr/local/bin/uv
ENV PATH="/root/.cargo/bin:${PATH}"

# Set up /io directory
RUN mkdir -p /io
RUN chown $NB_USER /io

USER $NB_USER

# Clone fidimag
ENV FIDIMAG_DIR=/io
WORKDIR /io
RUN git clone https://github.com/computationalmodelling/fidimag.git
WORKDIR /io/fidimag

# Install third party libraries from source
RUN bash bin/install-fftw.sh
RUN bash bin/install-sundials.sh

# Build fidimag using uv
RUN uv sync --frozen

ENV LD_LIBRARY_PATH=/io/fidimag/local/lib \
    OMP_NUM_THREADS=1

# https://github.com/conda-forge/matplotlib-feedstock/issues/36
RUN conda install --quiet --yes icu

# Verify installation
RUN uv run python -c "import fidimag; print(f'Fidimag {fidimag.__version__} installed')"

# /io will be mounted from the host system
WORKDIR /io

