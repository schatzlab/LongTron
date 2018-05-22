#!/bin/bash
#parameters:
#1: file of Target junctions
echo "lr_min_reads	target_min_reads	exact_matches	all_lr_jxs	all_target_jxs	lr_recall	target_recall" > compare.tsv
for i in 1 2 5 7 10 15 20 50 100 500 1000
do
	/bin/bash -x ./compare_junctions_target_filtered.sh NA12878.bam.perl.introns.merged.sorted $i $1 $i >> compare.tsv
done
