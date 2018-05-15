#!/bin/bash

#first filter lr jxs at given read count threshold
cat $1 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[3] < '${2}'); print "$f\n";' | cut -f 1,2,3 | sort -u > $1.filtered.${2}

#next filter Snaptron jxs at given total (sum over all samples) read count threshold
zcat $3 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[13] < '${4}'); print "$f\n";' | cut -f2,3,4 | sort -u > $3.filtered.${4}

comm -1 -2 $1.filtered.${2} $3.filtered.${4} | wc -l | perl -ne 'print "".'${2}'."\t".'${4}'."\t$_\n";' 
