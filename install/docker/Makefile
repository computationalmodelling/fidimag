# to build a new docker image
build:
	time docker build -t fangohr/fidimag:latest .

# to run new image
run:
	docker run -ti fangohr/fidimag bash
	# try 'ipython' and then 'import fidimag'

login:
	docker login

# to push the new docker image to dockerhub (need to login first)
push:
	docker push fangohr/fidimag:latest

# to fetch image to local machine
pull:
	docker pull fangohr/fidimag:latest
