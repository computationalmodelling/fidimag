# Docker

Run fidimag with Python 2 under Docker.

## Using the docker container

There is a fidimag container available under `fidimag/minimal-py2`.

To use it, you can try this:

1. Install docker

2. Pull the container onto your machine using

       docker pull fidimag/minimal-py2

3. Start the container using

       docker run -ti fidimag/minimal-py2

   This command should show a bash prompt inside the docker container:

    <pre>
    bin:docker fangohr$ docker run -ti fidimag/minimal-py2
    fidimag@38fdd2a0feb4:~$
    </pre>

One way to test the installation is to run the unit tests:

<pre>
fidimag@38fdd2a0feb4:~$ cd fidimag/tests/
fidimag@38fdd2a0feb4:~/fidimag/tests$ py.test -v 
========== test session starts =======================
platform linux2 -- Python 2.7.6 -- pytest-2.5.1 -- /usr/bin/python
collected 63 items 

field_test.py:7: test_initialise_scalar PASSED
field_test.py:13: test_initialise_vector PASSED
test_2dpbc_cube.py:10: test_compute_field PASSED
test_anis.py:7: test_anis PASSED
.
.
.
</pre>

## Shortcomings

- need to share directory between host and container (MOUNT or VOLUME)

## Creating the docker container

The `Makefile` in this directory shows the relevant targets (build, login, push)
to create a new container and push it to the the cloud.
