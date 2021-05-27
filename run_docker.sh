#!/bin/bash
docker run -v $(pwd):/work/jec-dnn -w /work/jec-dnn --user $(id -u):$(id -g) --gpus all -it jec-dnn