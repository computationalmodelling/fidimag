# to build a new docker image
build:
	time docker build -t fidimag/notebook:latest .

# to run new image
run: build
	docker run -v `pwd`:/io -d -p 30008:8888 fidimag/notebook

# to push the new docker image to dockerhub (need to login first)
push: build
	docker push fidimag/notebook:latest

# to fetch image to local machine
pull:
	docker pull fidimag/notebook:latest
