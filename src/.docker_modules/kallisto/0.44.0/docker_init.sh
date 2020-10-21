#!/bin/sh
docker pull lbmc/kallisto:0.44.0
docker build src/.docker_modules/kallisto/0.44.0 -t 'lbmc/kallisto:0.44.0'
docker push lbmc/kallisto:0.44.0
