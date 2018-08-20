#!/bin/bash
set -o pipefail -o nounset -o errexit 

#only pull out internal exons, therefore start/end/single exons are ignored
samtools view -F 2308 ${1} | cut -f 1,2,3,4,6 | { perl -ne 'BEGIN { $b=0; } chomp; ($tname,$flag,$c,$s,$f)=split(/\t/,$_); print STDERR "$b\t$tname\n"; my $o = (int($flag) & 0x10)?"-":"+"; $r=$s; $pr=$r; $first_exon=1; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N" && !$first_exon) { print "$c\t$pr\t".($r-1)."\t$b\t$o\n"; $pr=$r+$i; } $first_exon=undef; $r+=$i; } $b++;' 2>&1 1>&3 | sort -k2,2 > ${1}.t2ids.tsv; } 3>&1 | sort -k1,1 -k2,2n -k3,3n -k4,4 | perl -ne 'chomp; $f=$_; ($c,$s,$e,$sid,$o)=split(/\t/,$f); if($pc && $pc eq $c && $ps == $s && $pe == $e && $po eq $o) { $sids.=",$sid"; $nsids++; next; } elsif($pc) { print "$pc\t$ps\t$pe\t$po\t$nsids\t$sids\n"; } $nsids=1; $pc=$c; $ps=$s; $pe=$e; $po=$o; $sids=$sid; END { if($pc) { print "$pc\t$ps\t$pe\t$po\t$nsids\t$sids\n"; } }' | bgzip > ${2}

