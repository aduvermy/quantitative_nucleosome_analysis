#!/bin/sh
docker pull lbmc/bioawk:1.0
docker build src/.docker_modules/bioawk/1.0 -t 'lbmc/bioawk:1.0'
docker push lbmc/bioawk:1.0
