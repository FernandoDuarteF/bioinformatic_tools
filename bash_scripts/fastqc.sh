#!/usr/bin/env bash
READS=$1
OutDir=$2

CUR_PATH=$PWD/$OutDir

mkdir -p $CUR_PATH
cd $CUR_PATH

fastqc $READS

#unzip *.zip

exit