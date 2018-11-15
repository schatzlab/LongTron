#!/bin/env bash
set -o pipefail -o nounset -o errexit 
scripts=`perl -e '$f="'${0}'"; $f=~s/\/[^\/]+$/\//; print "$f\n";'`

#e.g. ../gencode.v28.basic.annotation.junctions4
junctions=$1
#e.g. ../gencode.v28.basic.annotation.transcripts_exon_count
transcript_exon_counts=$2
#if transcript names coming from the simulations need to have their suffixes trimed, pass a 1
fix_transcripts=$3

cd fl.20.all.0
for i in {0..4}; do cd ../fl.20.all.${i} && /bin/bash -x $scripts/junction_wiggle.sh trans_sim10.fl.junctions $junctions 20 trans_sim10.fl.bam.jxs.t2ids.tsv ; done

#now run the following
/bin/bash -x $scripts/compare_matching.sh $transcript_exon_counts $fix_transcripts > compare_matching.sh.run 2>&1
#then build features and use class2transcripts.tsv to map categories
