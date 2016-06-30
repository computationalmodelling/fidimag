# Docker

First steps towards running fidimag through Docker.

## Using the docker container

There is a fidimag container available under `fangohr/fidimag`.

To use it, you can try this:

1. Install docker

2. Pull the container onto your machine using

`docker pull fangohr/fidimag`

3. Start the container using

`docker run -ti fangohr/fidimag`

## Creating the docker container

## Shortcomings

- need to share directory between host and container (MOUNT or VOLUME)

- make use of jupyter notebook work (as we have seen in FEniCS course): connect port in container to port on host machine (EXPOSE)

## Creating the docker container

The `Makefile` in this directory shows the relevant targets (build, login, push) to create a new container and push it to the the cloud.





