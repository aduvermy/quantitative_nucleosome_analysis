#!/bin/sh
docker pull lbmc/deeptools:3.0.2
docker build src/.docker_modules/deeptools/3.0.2 -t 'lbmc/deeptools:3.0.2'
docker push lbmc/deeptools:3.0.2
