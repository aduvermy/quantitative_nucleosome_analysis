#!/bin/sh
docker pull lbmc/kallisto:0.43.1
docker build src/.docker_modules/kallisto/0.43.1 -t 'lbmc/kallisto:0.43.1'
docker push lbmc/kallisto:0.43.1
