#!/bin/bash

#e.g. ../gencode.v28.basic.annotation.junctions4
ANNOTATED_JXS=$1

for i in 0 1 2 3 4; do
	cd fl.20.all.${i}/
	/bin/bash -x ../../junction_wiggle.sh trans_sim10.fl.junctions $ANNOTATED_JXS 20 trans_sim10.fl.bam.jxs.t2ids.tsv
	cd ../
done
