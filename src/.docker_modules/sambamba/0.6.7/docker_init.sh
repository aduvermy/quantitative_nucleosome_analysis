#!/bin/sh
docker pull lbmc/sambamba:0.6.7
docker build src/.docker_modules/sambamba/0.6.7 -t 'lbmc/sambamba:0.6.7'
docker push lbmc/sambamba:0.6.7
