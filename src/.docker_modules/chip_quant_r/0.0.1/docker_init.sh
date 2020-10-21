#!/bin/sh
docker pull lbmc/chip_quant_r:0.0.1
docker build src/.docker_modules/chip_quant_r/0.0.1 -t 'lbmc/chip_quant_r:0.0.1'
docker push lbmc/chip_quant_r:0.0.1
