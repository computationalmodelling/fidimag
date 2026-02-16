FROM ubuntu:22.04

# Minimal Python 3 image for running fidimag
ARG DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get -y update
RUN apt-get -y install \
    build-essential cmake gcc g++ \
    libatlas-base-dev \
    python3 python3-dev \
    wget curl git ca-certificates

# Install uv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh && \
    ln -s /root/.cargo/bin/uv /usr/local/bin/uv
ENV PATH="/root/.cargo/bin:${PATH}"

# Set up fidimag user (do this before cloning to set ownership correctly)
RUN useradd -m -s /bin/bash -G sudo fidimag && \
    echo "fidimag:docker" | chpasswd && \
    echo "fidimag ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

# where to install source
ENV FIDIMAG_HOME=/home/fidimag

RUN mkdir -p $FIDIMAG_HOME
WORKDIR $FIDIMAG_HOME

# Clone fidimag
RUN git clone https://github.com/computationalmodelling/fidimag.git
WORKDIR $FIDIMAG_HOME/fidimag

# Install third party libraries from source
RUN bash bin/install-fftw.sh
RUN bash bin/install-sundials.sh

# Build fidimag using uv
RUN uv sync --frozen

ENV LD_LIBRARY_PATH=$FIDIMAG_HOME/fidimag/local/lib \
    OMP_NUM_THREADS=1

# Check that tests run okay
RUN uv run pytest tests/ -v -m "not slow and not run_oommf"

# Set up permissions
RUN chown -R fidimag $FIDIMAG_HOME

# For bind mounts from host
RUN mkdir /io
RUN chown -R fidimag /io

USER fidimag
RUN touch $FIDIMAG_HOME/.sudo_as_admin_successful

WORKDIR /io
CMD ["/bin/bash", "-i"]
