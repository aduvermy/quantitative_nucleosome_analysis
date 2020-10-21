#!/bin/sh
docker pull lbmc/rsem:1.3.0
docker build src/.docker_modules/rsem/1.3.0 -t 'lbmc/rsem:1.3.0'
docker push lbmc/rsem:1.3.0
