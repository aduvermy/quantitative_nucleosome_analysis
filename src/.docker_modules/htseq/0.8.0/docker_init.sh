#!/bin/sh
docker pull lbmc/htseq:0.8.0
docker build src/.docker_modules/htseq/0.8.0 -t 'lbmc/htseq:0.8.0'
docker push lbmc/htseq:0.8.0
