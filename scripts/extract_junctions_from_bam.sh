#!/bin/bash

#extract out properly mapped reads that are primary alignments, then parse out junctions in alignment,
#then merge all unique junction coordinates and record all supporting read indexes
samtools view -F 2308 $1 | cut -f 3,4,6 | perl -ne 'BEGIN { $b=0; } chomp; ($c,$s,$f)=split(/\t/,$_); $r=$s; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N") { print "$c\t$r\t".($r+$i-1)."\t$b\n"; } $r+=$i; } $b++;' | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; $f=$_; ($c,$s,$e,$sid)=split(/\t/,$f); if($pc && $pc eq $c && $ps == $s && $pe == $e) { $sids.=",$sid"; $nsids++; next; } elsif($pc) { print "$pc\t$ps\t$pe\t$nsids\t$sids\n"; } $nsids=1; $pc=$c; $ps=$s; $pe=$e; $sids=$sid; END { if($pc) { print "$pc\t$ps\t$pe\t$nsids\t$sids\n"; } }' > $1.introns
