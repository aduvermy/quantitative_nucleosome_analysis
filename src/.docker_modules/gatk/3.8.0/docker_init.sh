#!/bin/sh
docker pull lbmc/gatk:3.8.0
docker build src/.docker_modules/gatk/3.8.0 -t 'lbmc/gatk:3.8.0'
docker push lbmc/gatk:3.8.0
