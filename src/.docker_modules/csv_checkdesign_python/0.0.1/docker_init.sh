#!/bin/sh
docker pull lbmc/csv_checkdesign_python:0.0.1
docker build src/.docker_modules/csv_checkdesign_python/0.0.1 -t 'lbmc/csv_checkdesign_python:0.0.1'
docker push lbmc/csv_checkdesign_python:0.0.1
