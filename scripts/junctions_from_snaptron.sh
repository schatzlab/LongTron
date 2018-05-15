#!/bin/bash

#only grabs minimal set of fields (coords, strand, reduced annotation status (0/1 for whole and individual SS, number of samples, number of total reads)
cat /dev/stdin | cut -f 2-4,6,7,10,11,13,14 | perl -ne 'chomp; @f=split(/\t/,$_); $f[5]=1 if($f[5] ne "0"); $f[6]=1 if($f[6] ne "0"); print "".join("\t",@f)."\n";' > $1.cut

