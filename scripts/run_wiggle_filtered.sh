#!/bin/bash
#parameters:
#1: file of Target junctions
#2: path to directory of already run exact, filtered matches (so we can get total wc numbers for filtered sets)
source ./k8_env.sh
paf=./paftools.js
wiggle=10
cat ${1} | perl -ne 'chomp; ($c,$s,$e,$o,$a1,$a2,$a3,$ns,$nr)=split(/\t/,$_); $i++; print "$c\t\texon\t$s\t$e\t$nr\t?\t\ttranscript_id "tid_$i";\n";' > ${1}.gtf
#run raw comparison with wiggle, this results in list of matching introns with read counts
$paf junceval -l ${wiggle} -e -p ${1}.gtf <(samtools view -h NA12878.bam) 2> convert.stderr | cut -f 1,2,3,6 | sort | uniq -c > ${1}.gtf.compare

echo "lr_min_reads	target_min_reads	exact_matches	all_lr_jxs	all_target_jxs	lr_recall	target_recall" > compare.tsv
for i in 1 2 5 7 10 15 20 50 100 500 1000
do
	/bin/bash -x ./compare_wiggle_filtered.sh $i ../${1}.gtf.compare $i ${2} >> compare.tsv
done

