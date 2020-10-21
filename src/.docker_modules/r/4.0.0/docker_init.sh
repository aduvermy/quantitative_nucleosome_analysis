#!/bin/sh
docker pull lbmc/r-base:4.0.0
docker build src/.docker_modules/r/4.0.0 -t 'lbmc/r-base:4.0.0'
docker push lbmc/r-base:4.0.0
