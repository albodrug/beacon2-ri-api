# shut down running container
sudo docker-compose down
# create docker image of network
sudo docker network create my-app-network
# up the container
sudo docker-compose up --build
