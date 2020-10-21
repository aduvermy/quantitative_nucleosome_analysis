#!/bin/sh
docker pull lbmc/liftover:357
docker build src/.docker_modules/liftover/357/ -t 'lbmc/liftover:357'
docker push lbmc/liftover:357
