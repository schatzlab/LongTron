#!/bin/bash
set -o pipefail -o nounset -o errexit 
#extracts all splice junctions from each read which is:
#1) properly aligned
#2) and it not secondary or supplementary
#params:
#1: path to samtools
#2: path to BAM file
#3: path to splice output file
#4: path to isoform output file
#5: working dir (to create fifo in)
#6: temp directory for sorting

${1} view -F 2308 ${2} | cut -f 2,3,4,6 | perl -ne 'BEGIN { $b=0; } chomp; ($flag,$c,$s,$f)=split(/\t/,$_); my $o = (int($flag) & 0x10)?"-":"+"; $r=$s; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N") { print "$c\t$r\t".($r+$i-1)."\t$b\t$o\n"; } $r+=$i; } $b++;' | sort -T ${6} -k1,1 -k2,2n -k3,3n -k5,5 -k4,4n | perl -ne 'chomp; $f=$_; ($c,$s,$e,$sid,$o)=split(/\t/,$f); if($pc && $pc eq $c && $ps == $s && $pe == $e && $po eq $o) { $sids.=",$sid"; $nsids++; next; } elsif($pc) { print "$pc\t$ps\t$pe\t$po\t$nsids\t$sids\n"; } $nsids=1; $pc=$c; $ps=$s; $pe=$e; $po=$o; $sids=$sid; END { if($pc) { print "$pc\t$ps\t$pe\t$po\t$nsids\t$sids\n"; } }' 2> ${3}.err | bgzip > ${3}

${1} view -F 2308 ${2} | cut -f1,2,3,4,6 | extract_isoforms.pl 2> ${4}.err | gzip > ${4}
