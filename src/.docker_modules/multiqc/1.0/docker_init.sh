#!/bin/sh
docker pull lbmc/multiqc:1.0
docker build src/.docker_modules/multiqc/1.0 -t 'lbmc/multiqc:1.0'
docker push lbmc/multiqc:1.0
