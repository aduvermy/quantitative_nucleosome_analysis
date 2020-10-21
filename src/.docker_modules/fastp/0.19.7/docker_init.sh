#!/bin/sh
docker pull lbmc/fastp:0.19.7
docker build src/.docker_modules/fastp/0.19.7 -t 'lbmc/fastp:0.19.7'
docker push lbmc/fastp:0.19.7
