#!/usr/bin/env bash

BIGBED=/data/kent_tools/bedToBigBed
#CHROMSZ=/data/kent_tools/hg38.chrom.sizes
#CHROMSZ=./SRR1163655.sorted.bam.chr_sizes

#all.matches
suffix=$1
#fl or nonfl (usually fl)
base_type=$2
#for making BB files to display on UCSC GB
BAM_BED_FILE=$3
#need the specific set for chromosomes, e.g. SRR1163655.sorted.bam.chr_sizes 
CHROMSZ=$4

for f in problem-free non-recurrent recurrent novel; do
    #column to compare with (nonfl==12)
    starting_col=`perl -e '$t="'${base_type}'"; $starting_col=9; $starting_col=3 if($t eq "nonfl"); print $starting_col;'`
    col=`perl -e '%h=("problem-free"=>0,"non-recurrent"=>1,"recurrent"=>2,"novel"=>3); print "".($h{"'${f}'"}+'${starting_col}');'`
    cat ${f}.${base_type}.${suffix} | perl -ne 'BEGIN { $col='${col}'-1; $starting_col = '${starting_col}'-1; } chomp; $f=$_; @f=split(/\t/,$f); $max=0; $max_col=-1; for($i=$starting_col; $i<($starting_col+4); $i++) { if($f[$i] >= $max) { $max=$f[$i]; $max_col=$i; } } if($max_col == $col) { print "$f\n"; }' > ${f}.all.matches.joined
    cut -f 2 ${f}.all.matches.joined > ${f}.all.matches.joined.read_names
    export LC_COLLATE=C && fgrep -f ${f}.all.matches.joined.read_names $BAM_BED_FILE | sort -k1,1 -k2,2n -k3,3n > ${f}.all.matches.joined.reads.bed
    ${BIGBED} ${f}.all.matches.joined.reads.bed $CHROMSZ ${f}.all.matches.joined.reads.bed.bb
done
wc -l *.joined > all.joined.info
