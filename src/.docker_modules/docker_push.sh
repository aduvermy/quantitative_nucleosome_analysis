#!/bin/sh
fd "Dockerfile" src/.docker_modules | perl -pe 's|.*docker_modules/(.*)/(.*)/Dockerfile|\1:\2|g' | awk '{system("docker push lbmc/"$0)}'
