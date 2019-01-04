#!/bin/bash
set -o pipefail -o nounset -o errexit 

#extract out properly mapped reads that are primary alignments, then parse out junctions in alignment,
#then merge all unique junction coordinates and record all supporting read indexes
samtools view -F 2308 ${1} | cut -f 1,2,3,4,6 | { perl -ne 'BEGIN { $b=0; } chomp; ($tname,$flag,$c,$s,$f)=split(/\t/,$_); print STDERR "$b\t$tname\n"; my $o = (int($flag) & 0x10)?"-":"+"; $r=$s; $jid=0; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N") { $jid++; print "$c\t$r\t".($r+$i-1)."\t$b.$jid\t$o\n"; } $r+=$i; } $b++;' 2>&1 1>&3 | sort -T ${3} -k2,2 > ${1}.jxs.t2ids.tsv; } 3>&1 | sort -T ${3} -k1,1 -k2,2n -k3,3n -k5,5 -k4,4n | perl -ne 'chomp; $f=$_; ($c,$s,$e,$sid,$o)=split(/\t/,$f); if($pc && $pc eq $c && $ps == $s && $pe == $e && $po eq $o) { $sids.=",$sid"; $nsids++; next; } elsif($pc) { print "$pc\t$ps\t$pe\t$po\t$nsids\t$sids\n"; } $nsids=1; $pc=$c; $ps=$s; $pe=$e; $po=$o; $sids=$sid; END { if($pc) { print "$pc\t$ps\t$pe\t$po\t$nsids\t$sids\n"; } }' 2> ${2}.err > ${2}
