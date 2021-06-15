#!/bin/bash

print_usage() {
  printf "Usage: 
  -o out directory, -c config directory, 
  -d data directory, -g gpu id(s)
"
}

while getopts "o:c:d:g:" arg; do
  case $arg in
    o) out_dir=$OPTARG;;
    c) config_dir=$OPTARG;;
    d) data_dir=$OPTARG;;
    g) gpus=$OPTARG;;
    *) print_usage
       exit 1;;
  esac
done

mkdir -p $out_dir

for f in $config_dir/*
do
    filename=$(basename -- $f)
    name=${filename%.*}
    python train.py -i $data_dir -o $out_dir/$name -c $f --gpus $gpus
    python plot.py -i $out_dir/$name -o $out_dir/$name
done
