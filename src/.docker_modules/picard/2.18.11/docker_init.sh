#!/bin/sh
docker pull lbmc/picard:2.18.11
docker build src/.docker_modules/picard/2.18.11 -t 'lbmc/picard:2.18.11'
docker push lbmc/picard:2.18.11
