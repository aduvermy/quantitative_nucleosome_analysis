#!/bin/sh
docker pull lbmc/star:2.7.3a
docker build src/.docker_modules/star/2.7.3a/ -t 'lbmc/star:2.7.3a'
docker push lbmc/star:2.7.3a
