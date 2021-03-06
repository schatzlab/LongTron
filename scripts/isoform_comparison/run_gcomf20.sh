#!/usr/bin/env bash
set -ex -o pipefail

dir=$1
gtf1=$2
gtf2=$3

mkdir -p $dir
cd $dir
ln -fs ../${gtf1}
ln -fs ../${gtf2}
../gcom -a -T -r $gtf1 --fuzz-length 20 $gtf2 > run 2>&1
cd ../
