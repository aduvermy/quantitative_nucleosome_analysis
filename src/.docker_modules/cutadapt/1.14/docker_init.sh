#!/bin/sh
docker pull lbmc/cutadapt:1.14
docker build src/.docker_modules/cutadapt/1.14 -t 'lbmc/cutadapt:1.14'
docker push lbmc/cutadapt:1.14
