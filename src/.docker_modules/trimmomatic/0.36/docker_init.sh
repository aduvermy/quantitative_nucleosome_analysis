#!/bin/sh
docker pull lbmc/trimmomatic:0.36
docker build src/.docker_modules/trimmomatic/0.36 -t 'lbmc/trimmomatic:0.36'
docker push lbmc/trimmomatic:0.36
