#!/bin/bash
#parameters:
#1: file of LR junctions w/ # of reads supporting in col. 4
#2: min # of reads filter for LR
#3: file of Target junctions
#4: total # of junctions in Target

#first filter lr jxs at given read count threshold
cat $1 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[3] < '${2}'); print "$f\n";' | cut -f 1,2,3 | sort -u > $1.filtered.${2}
lr=`wc -l $1.filtered.${2} | cut -d" " -f 1`

comm -1 -2 $1.filtered.${2} $3 | wc -l | perl -ne 'chomp; $s=$_; $lr_per=$s/'${lr}'; $snap_per=$s/'${4}'; print "'${2}'\t1\t$s\t'${lr}'\t'${4}'\t"; printf("%.3f\t%.3f\n",$lr_per,$snap_per);' 
