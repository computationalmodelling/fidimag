# to build a new docker image
build:
	time docker build -t fidimag/notebook:latest .

# to run new image
run:
	docker run -ti fidimag/notebook bash

login:
    docker login

# to push the new docker image to dockerhub (need to login first)
push: build login
	docker push fidimag/notebook:latest

# to fetch image to local machine
pull:
	docker pull fidimag/notebook:latest
