# Docker

First steps towards running fidimag through Docker.

## Using the docker container

There is a fidimag container available under `fangohr/fidimag`.

To use it, you can try this:

1. Install docker

2. Pull the container onto your machine using

       docker pull fidimag/notebook

3. Start the container using

       docker run -d -p 30008:8888 fidimag/notebook

   This will start a notebook server. You can see it in your browser at
http://localhost:30008/ on Linux, or http://192.168.99.100:30008/ on Mac (you may
need to change this IP address if your docker VM is at a different address).

To run a shell instead of the notebook server, run:

    docker run -ti fidimag/notebook bash

## Shortcomings

- need to share directory between host and container (MOUNT or VOLUME)


## Creating the docker container

The `Makefile` in this directory shows the relevant targets (build, login, push)
to create a new container and push it to the the cloud.





