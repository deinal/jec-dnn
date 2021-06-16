#!/bin/bash

while getopts "d:" arg; do
  case $arg in
    d) dir=$OPTARG;;
    *) echo "Specify directory using -d";;
  esac
done

for subdir in $dir/*
do
  find $subdir -name *.png -mindepth 2 -delete
  python plot.py -i $subdir -o $subdir
done