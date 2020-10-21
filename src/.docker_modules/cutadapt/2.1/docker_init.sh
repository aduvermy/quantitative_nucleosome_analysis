#!/bin/sh
docker pull lbmc/cutadapt:2.1
docker build src/.docker_modules/cutadapt/2.1 -t 'lbmc/cutadapt:2.1'
docker push lbmc/cutadapt:2.1
