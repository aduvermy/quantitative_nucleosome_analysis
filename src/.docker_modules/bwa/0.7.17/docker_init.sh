#!/bin/sh
docker pull lbmc/bwa:0.7.17
docker build src/.docker_modules/bwa/0.7.17 -t 'lbmc/bwa:0.7.17'
docker push lbmc/bwa:0.7.17
