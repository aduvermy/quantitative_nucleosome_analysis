#!/bin/sh
docker pull lbmc/file_handle:0.1.1
docker build src/.docker_modules/file_handle/0.1.1 -t 'lbmc/file_handle:0.1.1'
docker push lbmc/file_handle:0.1.1
