.PHONY : docker_run docker_shell docker_rstudio docker_start \
docker_stop docker_rm_container docker_build docker_push docker_login \
docker_rm_image docker_pull

docker_run:
	docker run -it --name $(CONTAINER_NAME) -w /data -v $(shell pwd):/data -d -P $(IMAGE_NAME)

docker_shell:
	docker exec -it $(CONTAINER_NAME) bash

docker_rstudio:
	$(INTERNET_BROWSER) $(shell docker port $(CONTAINER_NAME) 8787) &

docker_start:
	docker start $(CONTAINER_NAME)

docker_stop:
	docker stop $(CONTAINER_NAME)

docker_rm_container:
	docker rm $(CONTAINER_NAME)

docker_build: Dockerfile
	docker build -t $(IMAGE_NAME) .

docker_push:
	docker push $(IMAGE_NAME)

docker_login:
	docker login --username=cayek --email=kevin-ca@club-internet.fr

docker_rm_image:
	docker rmi -f $(IMAGE_NAME)

docker_pull:
	docker pull $(IMAGE_NAME)
