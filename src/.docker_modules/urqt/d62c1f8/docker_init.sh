#!/bin/sh
docker pull lbmc/urqt:d62c1f8
docker build src/.docker_modules/urqt/d62c1f8 -t 'lbmc/urqt:d62c1f8'
docker push lbmc/urqt:d62c1f8
