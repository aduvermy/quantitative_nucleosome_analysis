#!/bin/sh
docker pull lbmc/fastqc:0.11.5
docker build src/.docker_modules/fastqc/0.11.5 -t 'lbmc/fastqc:0.11.5'
docker push lbmc/fastqc:0.11.5
