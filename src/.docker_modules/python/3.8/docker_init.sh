#!/bin/sh
docker pull lbmc/python:3.8
docker build src/.docker_modules/python/3.8 -t 'lbmc/python:3.8'
docker push lbmc/python:3.8
