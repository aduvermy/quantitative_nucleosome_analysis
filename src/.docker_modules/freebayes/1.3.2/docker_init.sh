#!/bin/sh
docker pull lbmc/freebayes:1.3.2
docker build src/.docker_modules/freebayes/1.3.2/ -t 'lbmc/freebayes:1.3.2'
docker push lbmc/freebayes:1.3.2
