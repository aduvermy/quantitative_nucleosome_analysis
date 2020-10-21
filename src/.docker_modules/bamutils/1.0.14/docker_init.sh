#!/bin/sh
docker pull lbmc/bamutils:1.0.14
docker build src/.docker_modules/bamutils/1.0.14 -t 'lbmc/bamutils:1.0.14'
docker push lbmc/bamutils:1.0.14
