# Docker

Run fidimag with Python 3 under Docker.

## Using the docker container

There is a fidimag container available under `fidimag/minimal-py3`.

To use it, you can try this:

1. Install docker

2. Pull the container onto your machine using

       docker pull fidimag/minimal-py3

3. Start the container using

       docker run -ti fidimag/minimal-py3

   This command should show a bash prompt inside the docker container:

    <pre>
    bin:docker fangohr$ docker run -v `pwd`:/io -ti fidimag/minimal-py3
    fidimag@38fdd2a0feb4:/io$
    </pre>

## Creating the docker container

The `Makefile` in this directory shows the relevant targets (build, login, push)
to create a new container and push it to the the cloud.
