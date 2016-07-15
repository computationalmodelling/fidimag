# to build a new docker image
build:
	time docker build -t fidimag/minimal-py2:latest .

# to run new image
run: build
	docker run -v `pwd`:/io -ti fidimag/minimal-py2
	# try 'ipython' and then 'import fidimag'

# to push the new docker image to dockerhub (need to login first)
push: build
	docker push fidimag/minimal-py2:latest

# to fetch image to local machine
pull:
	docker pull fidimag/minimal-py2:latest
