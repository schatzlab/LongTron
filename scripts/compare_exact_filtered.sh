#!/bin/bash
#parameters:
#1: file of LR junctions w/ # of reads supporting in col. 4
#2: min # of reads filter for LR
#3: file of Target junctions w/ # of reads supporting in col. 9 (either for a single sample or all samples if from Snaptron)
#4: min # of reads filter for Target

#first filter lr jxs at given read count threshold
#cat $1 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[4] < '${2}'); print "$f\n";' | cut -f 1,2,3 | sort -u > lr.filtered.${2}
ln -s ../lr_filtered/lr.filtered.${2}
lr=`wc -l lr.filtered.${2} | cut -d" " -f 1`

#next filter pre-formatted target jxs at given total (sum over all samples) read count threshold
cat $3 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[8] < '${4}'); print "$f\n";' | cut -f1,2,3 | sort -u > target.filtered.${4}
target=`wc -l target.filtered.${4} | cut -d" " -f 1`


comm -1 -2 lr.filtered.${2} target.filtered.${4} | wc -l | perl -ne 'chomp; $s=$_; $lr_per=$s/'${lr}'; $target_per=$s/'${target}'; print "'${2}'\t'${4}'\t$s\t'${lr}'\t'${target}'\t"; printf("%.3f\t%.3f\n",$lr_per,$target_per);' 
