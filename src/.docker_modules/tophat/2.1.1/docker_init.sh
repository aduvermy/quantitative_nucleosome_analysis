#!/bin/sh
docker pull lbmc/tophat:2.1.1
docker build src/.docker_modules/tophat/2.1.1 -t 'lbmc/tophat:2.1.1'
docker push lbmc/tophat:2.1.1
