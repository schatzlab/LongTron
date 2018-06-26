#!/bin/bash
BEDTOOLS='/data/bedtools2/bin/bedtools'

WIGGLE=${1}
cat pre_novel_exons | perl -ne 'BEGIN { $w='${WIGGLE}'; } chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$exons,$len)=split(/\t/,$f); $s1=$s-($w+1); $s2=$s+$w; $e1=$e-($w+1); $e2=$e+$w; print "".join("\t",($c,$s1,$s2,$o,$nr,$reads,$pac,$nlr,$dist,"$s-$e",0,$exons,$len))."\n"; print "".join("\t",($c,$e1,$e2,$o,$nr,$reads,$pac,$nlr,$dist,"$s-$e",1,$exons,$len))."\n";' | sort -k1,1 -k2,2n -k3,3n > pre_novel_exons.split_ends.${WIGGLE}.bed

#find non-overlaps with list of all short read exons looking for ends separately (to allow for distant containments/overlaps)
$BEDTOOLS intersect -sorted -a pre_novel_exons.split_ends.${WIGGLE}.bed -b <(zcat ../gtex_sra_junctions.split.bed.bgz) -v > novel_exons_bt.w${WIGGLE}.raw

cat novel_exons_bt.w${WIGGLE}.raw | perl -ne 'BEGIN { %h; } chomp; $f=$_; ($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$coord,$t,$exons,$len)=split(/\t/,$f); $tnot=!$t; $k=$coord."-".$tnot; if($h{$k}) { print "".$h{$k}."\n"; delete($h{$k}); next; } $k=$coord."-".$t; $coord=~/^(\d+)-(\d+)$/; $c1=$1; $c2=$2; $h{$k}=join("\t",($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$exons,$len));' | egrep -v -e '(_alt)|(_random)|(_decoy)' > novel_exons_bt.w${WIGGLE}.both_ends

cat novel_exons_bt.w${WIGGLE}.raw | perl -ne 'BEGIN { %h; } chomp; $f=$_; ($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$coord,$t,$exons,$len)=split(/\t/,$f); $k=$coord; $coord=~/^(\d+)-(\d+)$/; $c1=$1; $c2=$2; if(!$h{$k}) { print "".join("\t",($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$exons,$len))."\n"; } $h{$k}=1;' | egrep -v -e '(_alt)|(_random)|(_decoy)' > novel_exons_bt.w${WIGGLE}.either_end

cat novel_exons_bt.w${WIGGLE}.either_end | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[8] == 0); if($f[6]==1) { print "$f\n"; } else { print STDERR "$f\n";}' > novel_exons_bt.w${WIGGLE}.either_end.pb 2> novel_exons_bt.w${WIGGLE}.either_end.npb


sort -k5,5nr -k9,9nr novel_exons_bt.w${WIGGLE}.either_end.pb > novel_exons_bt.w${WIGGLE}.either_end.pb.sorted
