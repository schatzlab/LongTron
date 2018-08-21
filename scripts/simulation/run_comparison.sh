#!/bin/bash

fltype=$1
wiggle=$2
cutoff=$3

mkdir ${fltype}.${wiggle}.${cutoff}
cd ${fltype}.${wiggle}.${cutoff}
ln -fs ../counts.sh ../exon_wiggle.sh ../gencode.v28.basic.annotation.exons.sorted ../${fltype}.${cutoff}.exons.clean ./
/bin/bash -x exon_wiggle.sh ${fltype}.${cutoff}.exons gencode.v28.basic.annotation.exons.sorted $wiggle

