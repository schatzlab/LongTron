#!/bin/bash
#creates the file which contains all matching (within wiggle)
#junctions between the source long read junctions and a target (e.g. gtex) junctions
#input:
#1: source long read junctions (sorted by sort -k1,1 -k2,2n -k3,3n): chr start end num_supporting_reads (e.g. NA12878.bam.perl.introns.merged.sorted)
#2: source target junctions (sorted by sort -k1,1 -k2,2n -k3,3n): chr start end name num_supporting_reads (e.g. gtex_junctions.bed)
WIGGLE=10
BEDTOOLS=~/bedtools2/bin/bedtools

$BEDTOOLS intersect -sorted -wo -a $1 -b $2 | perl -ne 'BEGIN { $MAX_DIFF='${WIGGLE}'; } chomp; $f=$_; ($c1,$s1,$e1,$nr1,$c2,$s2,$e2,$n2,$nr2)=split(/\t/,$f); $d1=abs($s1-$s2); $d2=abs($e1-$e2); if($d1 <= $MAX_DIFF && $d2 <= $MAX_DIFF) { print "$f\n"; }' | sort -k1,1 -k2,2n -k3,3n -k9,9nr | perl filter_for_best_single_wiggle_match.pl > wiggle.lr.target.intersect.single 
