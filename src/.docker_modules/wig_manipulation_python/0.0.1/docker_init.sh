#!/bin/sh
docker pull lbmc/wig_manipulation_python:0.0.1
docker build src/.docker_modules/wig_manipulation_python/0.0.1 -t 'lbmc/wig_manipulation_python:0.0.1'
docker push lbmc/wig_manipulation_python:0.0.1
