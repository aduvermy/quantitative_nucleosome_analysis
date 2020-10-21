#!/bin/sh
docker pull lbmc/salmon:0.8.2
docker build src/.docker_modules/salmon/0.8.2 -t 'lbmc/salmon:0.8.2'
docker push lbmc/salmon:0.8.2
