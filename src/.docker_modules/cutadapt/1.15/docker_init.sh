#!/bin/sh
docker pull lbmc/cutadapt:1.15
docker build src/.docker_modules/cutadapt/1.15 -t 'lbmc/cutadapt:1.15'
docker push lbmc/cutadapt:1.15
