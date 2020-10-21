#!/bin/sh
docker pull lbmc/macs2:2.1.2
docker build src/.docker_modules/macs2/2.1.2 -t 'lbmc/macs2:2.1.2'
docker push lbmc/macs2:2.1.2
