#!/bin/sh

## Small script for concatenating read.bed files in each inversion direcotory, outputs them into main inversion dir as 'cat_reads.bed.gz'
## This way you can upload the one file to the UCSC genome browser instead of uploading all of them one at a time (killer I know)


for file in HOTSPOT_EVENTS/*/
do
  #echo $file
  dir=$file
  base="reads/"
  cat=$dir$base
  #echo $cat/*
  #rm $dir/cat_reads.bed
  cat $cat/* > $dir/cat_reads.bed.gz
  rm $cat/*
  
done
