#!/bin/sh
docker pull lbmc/r:3.5.3
docker build src/.docker_modules/r/3.5.3 -t 'lbmc/r:3.5.3'
docker push lbmc/r:3.5.3
