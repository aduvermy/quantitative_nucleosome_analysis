#!/bin/sh
docker pull lbmc/gatk:4.0.8.1
docker build src/.docker_modules/gatk/4.0.8.1 -t 'lbmc/gatk:4.0.8.1'
docker push lbmc/gatk:4.0.8.1
