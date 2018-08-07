#!/bin/bash

CUFFCOMP=~/cufflinks/cuffcompare
#ANNOTATION=gencode.v28.basic.annotation.gtf

BAM=$1
ANNOTATION=$2

/data3/LongReadRNA/scripts/pipeline/extract_exons.sh ${BAM} ${BAM}.exons
samtools view -F 2308 ${BAM} | cut -f1,2,3,4,6 | /data3/LongReadRNA/scripts/isoforms/extract_isoforms.pl 2> ${BAM}.isoforms.err > ${BAM}.cuff.gtf
$CUFFCOMP -C -F -V ${BAM}.cuff.gtf -o ${BAM}.cuff.collapsed.gtf > ${BAM}.cuff.collapse.processed_run 2>&1
$CUFFCOMP -r $ANNOTATIOn -C -F -V ${BAM}.cuff.collapsed.gtf.combined.gtf -o ${BAM}.combined_cuff.vs.orig_annotation > ${BAM}.combined_cuff.vs.orig_annotation.run 2>&1
#to be used in creating the BigBed file for visualization on UCSC GB
cat ${BAM}.cuff.collapsed.gtf.combined.gtf | perl -ne 'chomp; $s=$_; $s=~/oId "([^\.]+)\./; $d=$1; print "$d\n";' | sort -u > ${BAM}.cuff.collapsed.gtf.combined.gtf.rids
/bin/bash -x make_bb_track.sh $BAM ${BAM}.cuff.collapsed.gtf.combined.gtf.rids 
