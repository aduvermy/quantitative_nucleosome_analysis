#!/bin/sh
docker pull lbmc/kallistobustools:0.24.4
docker build src/.docker_modules/kallistobustools/0.24.4 -t 'lbmc/kallistobustools:0.24.4'
docker push lbmc/kallistobustools:0.24.4
