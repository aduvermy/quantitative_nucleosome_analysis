#!/bin/sh
docker pull lbmc/samblaster:0.1.24
docker build src/.docker_modules/samblaster/0.1.24 -t 'lbmc/samblaster:0.1.24'
docker push lbmc/samblaster:0.1.24
