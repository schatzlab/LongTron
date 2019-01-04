#!/bin/bash
set -o pipefail -o nounset -o errexit 
#extracts all splice junctions from each read which is:
#1) properly aligned
#2) and it not secondary or supplementary
#params:
#1: path to BAM file
#2: path to intron/exon output file prefix
#3: path to isoform output file
#4: path to directory to use as temporary storage for sorting

scripts=`perl -e '@f=split(/\//,"'${0}'"); pop(@f); print "".join("/",@f)."\n";'`

#extract introns
/bin/bash -x ${scripts}/extract_junctions.sh ${1} ${2}.junctions ${4}

#extract exons
/bin/bash -x ${scripts}/extract_exons.sh ${1} ${2}.exons ${4}

samtools view -F 2308 ${1} | cut -f1,2,3,4,6 | ${scripts}/extract_isoforms.pl 2> ${3}.err > ${3}
