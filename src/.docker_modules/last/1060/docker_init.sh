#!/bin/sh
docker pull lbmc/last:1060
docker build src/.docker_modules/last/1060/ -t 'lbmc/last:1060'
docker push lbmc/last:1060
