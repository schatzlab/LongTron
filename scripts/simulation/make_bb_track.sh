#!/bin/bash

BIGBED=/data/kent_tools/bedToBigBed
CHROMSZ=/data/kent_tools/hg38.chrom.sizes

#read indices to pullout
BAM=$1
RFILE=$2
samtools view -F 2308 $BAM | cut -f 2,3,4,6 | perl -ne 'BEGIN { $b=0; open(IN,"<'${RFILE}'"); while($line=<IN>) { chomp($line); $h{$line}=1;} close(IN); } chomp; if(!$h{$b}) { $b++; next; } ($flag,$c,$s,$f)=split(/\t/,$_); my $o = (int($flag) & 0x10)?"-":"+"; $r=$s; $pr=$r; print "$c\t".($pr-1); $line=""; $bst=""; $bsz=""; $bc=1; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N") { $bc++; $bst.=($pr-$s).","; $bsz.=($r-$pr).","; $pr=$r+$i; } $r+=$i; } $bst.=($pr-$s); $bsz.=($r-$pr); print "\t".($r-1)."\tr$b\t0\t$o\t".($s-1)."\t".($r-1)."\t255,0,0\t$bc\t$bsz\t$bst\n"; $b++;' | sort -k1,1 -k2,2n -k3,3n > ${BAM}.bed
${BIGBED}  ${BAM}.bed $CHROMSZ ${BAM}.bed.bb
