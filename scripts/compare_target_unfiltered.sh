#!/bin/bash
#parameters:
#1: file of Target junctions
target_f=$1
target_wc=`wc -l $target_f | cut -d" " -f 1`
echo "lr_min_reads	snap_min_reads	exact_matches	all_lr_jxs	all_gtex_jxs	lr_recall	target_recall" > compare.tsv
for i in 1 10 20 50 100 500 1000
do
	/bin/bash -x ./compare_junctions_target_unfiltered.sh NA12878.bam.perl.introns.merged.sorted $i $target_f $target_wc >> compare.tsv
done
