#!/bin/sh
docker pull lbmc/umi_tools:1.0.0
docker build src/.docker_modules/umi_tools/1.0.0/ -t 'lbmc/umi_tools:1.0.0'
docker push lbmc/umi_tools:1.0.0
