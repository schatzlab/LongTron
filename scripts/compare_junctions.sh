#!/bin/bash

#first filter lr jxs at given read count threshold
cat $1 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[3] < '${2}'); print "$f\n";' | cut -f 1,2,3 | sort -u > $1.filtered.${2}
lr=`wc -l $1.filtered.${2} | cut -d" " -f 1`

#next filter pre-cut Snaptron jxs at given total (sum over all samples) read count threshold
cat $3 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[8] < '${4}'); print "$f\n";' | cut -f1,2,3 | sort -u > $3.filtered.${4}
snap=`wc -l $3.filtered.${4} | cut -d" " -f 1`

comm -1 -2 $1.filtered.${2} $3.filtered.${4} | wc -l | perl -ne 'chomp; $s=$_; $lr_per=$s/'${lr}'; $snap_per=$s/'${snap}'; print "'${2}'\t'${4}'\t$s\t'${lr}'\t'${snap}'\t"; printf("%.3f\t%.3f\n",$lr_per,$snap_per);' 
