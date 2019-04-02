#!/usr/bin/env bash
set -o pipefail -o nounset -o errexit 

dir=$(dirname $0)
annotation=/data7/schatz/comparisons/G029.gtf.introns.filtered_sorted.bed.formatted
snaptron=/data7/schatz/comparisons/srav2_gtex_tcga.junctions.sorted.bed.bgz
#e.g. NA12878-DirectRNA.bam.junctions.clean.sorted.bed
jxs=$1
#e.g. "oxford"
tech=$2

#cat $annotation | perl ${dir}/jx_exact.pl $jxs > ${jxs}.annotation_exact.tsv 2> ${jxs}.annotation_exact.tids.tsv
#zcat $snaptron | perl ${dirname}/jx_exact.pl $jxs > ${jxs}.snaptron_exact.tsv 2> ${jxs}.snaptron_exact.tids.tsv
#/usr/bin/time -v bedtools intersect -sorted -wao -a NA12878-DirectRNA.bam.junctions.clean.sorted.bed -b <(cat schatz/comparisons/G029.gtf.introns.filtered_sorted.bed.formatted) | LongReadRNA/scripts/simulation/classification/pbt -l 11 > oxford_vs_g029.wiggle.out 2>oxford_vs_g029.wiggle.out.err

type="annotation"
for f in $annotation $snaptron ; do
    #need total for finding percentage wiggle
    total=`cat $jxs | wc -l`
    zcat $f | perl ${dir}/jx_exact.pl $jxs > ${jxs}.${type}_exact.tsv 2> ${jxs}.${type}_exact.tids.tsv
    #get stats from exact matching
    exact=`head -1 ${jxs}.${type}_exact.tsv`
    bedtools intersect -sorted -wao -a $jxs -b <(zcat $f) | ${dir}/pbt -l 11 > ${tech}_vs_${type}.wiggle.out 2> ${tech}_vs_${type}.wiggle.out
    #find how many are annotated/snaptron in wiggle output
    wiggle=`cat ${tech}_vs_${type}.wiggle.out | egrep -e '	1$' | wc -l`
    perl -e '$t="'$total'"; $e="'$exact'"; $w="'$wiggle'"; $type="'$type'"; $tech="'$tech'"; $wp=($wigggle/$total); printf("$tech\t$type\t$exact\t$wiggle (%.0f\%)\n",100*$wp);'
    type="snaptron"
done
