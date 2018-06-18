#!/bin/bash
#1: raw read based isoforms sorted (e.g. NA12878-DirectRNA.raw_isoforms.tsv.sorted.bgz)

zcat $1 | egrep -e '	\+	' | sort -k2,2 -k6,6n -k4,4n > ${1}.plus.s
zcat $1 | egrep -e '	-	' | sort -k2,2 -k4,4n -k6,6n > ${1}.minus.s

#clustering is currently only by 1) chrom, 2) strand and 3) shared end
#further processing is needed to split clusters where junctions aren't near each other
#also, this throws out all non-spliced (singleton) reads
python cluster_and_collapse_isoforms.py ${1}.plus.s > ${1}.plus.s.clustered
python cluster_and_collapse_isoforms.py ${1}.minus.s > ${1}.minus.s.clustered
