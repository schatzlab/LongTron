#!/bin/bash
#parameters:

#1: path to directory of already run wiggle base intersection between lr and target
#(this is output from create_wiggle_base_intersection.sh)

#2: path to directory of already-run exact-matching-filtered version of this target
#this is for calculting the wc of the filtered version at each level for the target

echo "lr_min_reads	target_min_reads	exact_matches	all_lr_jxs	all_target_jxs	lr_recall	target_recall" > compare.tsv
for i in 1 2 5 7 10 15 20 50 100 500 1000
do
	/bin/bash -x ./compare_wiggle_filtered.sh $i ${1} $i ${2} >> compare.tsv
done

