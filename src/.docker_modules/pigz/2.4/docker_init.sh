#!/bin/sh
docker pull lbmc/pigz:2.4
docker build src/.docker_modules/pigz/2.4 -t 'lbmc/pigz:2.4'
docker push lbmc/pigz:2.4
