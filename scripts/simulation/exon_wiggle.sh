#!/bin/bash

BT=/data/bedtools2/bin/bedtools

EXONS=$1
#refseq_exons1.nostrand.split_ends.bed
#gencode_refseq_exons.nostrand.sorted
#gencode.v28.basic.annotation.exons.sorted
ANNOT_EXONS=$2
WIGGLE=$3

#filter out any single read supported exon
#egrep -v -e '(HLA)|(alt)|(random)|(decoy)|(chrUn)' ${exons} > ${EXONS}.clean

#totally novel with wiggle (no split ends complications)
$BT window -a <(cat ${EXONS}.clean | perl -ne 'chomp; ($c,$s,$e,$o,$nr,$reads)=split(/\t/,$_); $s--; print "".join("\t",($c,$s,$e,$o,$nr,$reads))."\n";') -b <( cat $ANNOT_EXONS | perl -ne 'chomp; ($c,$s,$e)=split(/\t/,$_); $s--; print "$c\t$s\t$e\n";') -w $WIGGLE -v | gzip > sim10.exons.clean.w${WIGGLE}.novel.gz


#split annotated exons to be both ends (but separate)
#cat $ANNOT_EXONS | perl -ne 'chomp; $f=$_; ($c,$s,$e)=split(/\t/,$f); $s1=$s-($w+1); $s2=$s+$w; $e1=$e-($w+1); $e2=$e+$w; print "".join("\t",($c,$s1,$s2))."\n"; print "".join("\t",($c,$e1,$e2))."\n";' | sort -k1,1 -k2,2n -k3,3n > $ANNOT_EXONS.split_ends.bed
#w/o window
cat $ANNOT_EXONS | perl -ne 'chomp; $f=$_; ($c,$s,$e)=split(/\t/,$f); $s1=$s-1; $s2=$s; $e1=$e-1; $e2=$e; print "".join("\t",($c,$s1,$s2))."\n"; print "".join("\t",($c,$e1,$e2))."\n";' | sort -k1,1 -k2,2n -k3,3n > $ANNOT_EXONS.split_ends.bed

#split aligned exons to be both ends (but separate)
#cat ${EXONS}.clean | perl -ne 'BEGIN { $w='${WIGGLE}'; } chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads)=split(/\t/,$f); next if($nr < 2); $s1=$s-($w+1); $s2=$s+$w; $e1=$e-($w+1); $e2=$e+$w; print "".join("\t",($c,$s1,$s2,$o,$nr,$reads,"$s-$e"))."\n"; print "".join("\t",($c,$e1,$e2,$o,$nr,$reads,"$s-$e"))."\n";' | sort -k1,1 -k2,2n -k3,3n > sim10.exons.split_ends.clean.nr2.${WIGGLE}.bed
#w/o window
cat ${EXONS}.clean | perl -ne 'BEGIN { $w='${WIGGLE}'; } chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads)=split(/\t/,$f); $s1=$s-1; $s2=$s; $e1=$e-1; $e2=$e; print "".join("\t",($c,$s1,$s2,$o,$nr,$reads,"$s-$e"))."\n"; print "".join("\t",($c,$e1,$e2,$o,$nr,$reads,"$s-$e"))."\n";' | sort -k1,1 -k2,2n -k3,3n > sim10.exons.split_ends.bed


#now check via ends

#use split ends to check either side (more flexible)
#$BT intersect -sorted -a sim10.exons.split_ends.clean.nr2.${WIGGLE}.bed -b $ANNOT_EXONS.split_ends.bed -v | gzip > sim10.refseq.exons.split_ends.clean.nr2.w${WIGGLE}.raw.gz
$BT window -a sim10.exons.split_ends.bed -b $ANNOT_EXONS.split_ends.bed -w $WIGGLE -v | gzip > sim10.exons.split_ends.bed.w${WIGGLE}.gz

#but now require both ends have to be novel
zcat sim10.exons.split_ends.bed.w${WIGGLE}.gz | perl -ne 'BEGIN { %h; } chomp; $f=$_; ($c,$c1,$c2,$o,$nr,$reads,$coord)=split(/\t/,$f); $k=$coord; if($h{$k}) { print "".$h{$k}."\n"; delete($h{$k}); next; } $k=$coord; $coord=~/^(\d+)-(\d+)$/; $c1=$1; $c2=$2; $h{$k}=join("\t",($c,$c1,$c2,$o,$nr,$reads));' | sort -k5,5nr -k1,1 -k2,2n -k3,3n > sim10.exons.split_ends.bed.w${WIGGLE}.both_ends

#get truely novel exons (on both ends), different from the step above, as that could allow containment where neither end is close to an existing exon's end
#whereas this requires 0 overlap (not containment)
#$BT intersect -sorted -a sim10.exons.split_ends.clean.nr2.${WIGGLE}.bed -b $ANNOT_EXONS -v | gzip > sim10.refseq.exons.split_ends.clean.nr2.w${WIGGLE}.raw.novel.gz
#$BT window -a sim10.exons.split_ends.clean.nr2.${WIGGLE}.bed -b $ANNOT_EXONS -w $WIGGLE -v | gzip > sim10.refseq.exons.split_ends.clean.nr2.w${WIGGLE}.raw.novel.gz

#either end can be novel, but print the end that we know *is* novel; also print the original (non-split, non-windowed) BED-formatted coordinates as primary, important for the rest of the steps
zcat sim10.exons.split_ends.bed.w${WIGGLE}.gz | perl -ne 'BEGIN { %h; } chomp; $f=$_; ($c,$c1,$c2,$o,$nr,$reads,$coord)=split(/\t/,$f); $k=$coord; $coord=~/^(\d+)-(\d+)$/; $c3=$c1; $c1=$1; $c2=$2; $d1=abs($c1-$c3); $d2=abs($c2-$c3); $c3=($d1>$d2?$c2:$c1); if(!$h{$k}) { print "".join("\t",($c,--$c1,$c2,$o,$nr,$reads,$c3))."\n"; } $h{$k}=1;' | sort -k5,5nr -k1,1 -k2,2n -k3,3n > sim10.exons.split_ends.bed.w${WIGGLE}.either_end

SPLIT=sim10.exons.split_ends.bed.w${WIGGLE}.either_end.sorted
sort -k1,1 -k2,2n -k3,3n sim10.exons.split_ends.bed.w${WIGGLE}.either_end > $SPLIT

#now separate into further categories those which have some form of overlap:

$BT window -a $SPLIT -b $ANNOT_EXONS -w $WIGGLE | uniq > ${SPLIT}.raw

#1) OVERLAPPING
#use the original full coord range (not split ends) for both of the sets of intervals this time to remove contained/overlapping exons
#$BT window -a $SPLIT -b $ANNOT_EXONS -w $WIGGLE | uniq | perl -ne 'BEGIN {$w='${WIGGLE}';} chomp; $f=$_; @f=split(/\t/,$f); ($s,$e)=($f[1],$f[2]); ($s2,$e2)=($f[8],$f[9]); $d1=abs($s2-$s); $d2=abs($e2-$e); $k="$c:$s-$e"; if($pk && $k ne $pk && !$bad{$pk} && $h{$pk}) { print $h{$pk}; } delete $h{$pk}; $pk=$k; if(!(($s < $s2 && $d1 > $w) && ($e > $e2 && $d2 > $w)) && ($d1 > $w || $d2 > $w) && !($s <= $s2 && $e >= $e2)) { $f=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/$1:$2-$3\t/;  $h{$k}="$f\t$d1\t$d2\n"; } else { $bad{$k}=1; } END { if($pk && !$bad{$pk} && $h{$pk}) { print $h{$pk}; } }' | sort -k3,3nr > w${WIGGLE}.refseq_gencode.overlapping_not_contained
cat ${SPLIT}.raw | perl -ne 'BEGIN {$w='${WIGGLE}';} chomp; $f=$_; @f=split(/\t/,$f); ($s,$e)=($f[1],$f[2]); ($s2,$e2)=($f[8],$f[9]); $d1=abs($s2-$s); $d2=abs($e2-$e); $k="$c:$s-$e"; if($pk && $k ne $pk && $h{$pk}) { print $h{$pk}; delete $h{$pk}; } $pk=$k; if( ($s < $s2 && $d1 > $w  && !($e > $e2 && $d2 > $w)) || ($e > $e2 && $d2 > $w && !($s < $s2 && $d1 > $w)) ) { $f=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/$1:$2-$3\t/;  $h{$k}="$f\t$d1\t$d2\n"; } END { if($pk && $h{$pk}) { print $h{$pk}; } }' | sort -k3,3nr > w${WIGGLE}.refseq_gencode.overlapping_not_contained

#2) CONTAINED (new exons within annotated)
#get the contained ones
#$BT window -a <(sort -k1,1 -k2,2n -k3,3n sim10.refseq.exons.split_ends.clean.nr2.w${WIGGLE}.either_end) -b $ANNOT_EXONS -w $WIGGLE | uniq | perl -ne 'BEGIN {$w='${WIGGLE}';} chomp; $f=$_; @f=split(/\t/,$f); ($s,$e)=($f[1],$f[2]); ($s2,$e2)=($f[8],$f[9]); $d1=abs($s2-$s); $d2=abs($e2-$e); $k="$c:$s-$e"; if($pk && $k ne $pk && !$bad{$pk} && $h{$pk}) { print $h{$pk}; } delete $h{$pk}; $pk=$k; if(!(!(($s < $s2 && $d1 > $w) || ($e > $e2 && $d2 > $w)) && ($d1 > $w || $d2 > $w) && !($s <= $s2 && $e >= $e2))) { $f=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/$1:$2-$3\t/;  $h{$k}="$f\t$d1\t$d2\n"; } else { $bad{$k}=1; } END { if($pk && !$bad{$pk} && $h{$pk}) { print $h{$pk}; } }' | sort -k3,3nr > w${WIGGLE}.refseq_gencode.contained_not_overlapping
cat ${SPLIT}.raw | perl -ne 'BEGIN {$w='${WIGGLE}';} chomp; $f=$_; @f=split(/\t/,$f); ($s,$e)=($f[1],$f[2]); ($s2,$e2)=($f[8],$f[9]); $d1=abs($s2-$s); $d2=abs($e2-$e); $k="$c:$s-$e"; if($pk && $k ne $pk && $h{$pk}) { print $h{$pk}; delete $h{$pk}; } $pk=$k; if( ($s >= $s2 && $e <= $e2) || ($s < $s2 && $d1 <= $w && !($e > $e2 && $d2 > $w)) || ($e > $e2 && $d2 <= $w && !($s < $s2 && $d1 > $w)) ) { $f=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/$1:$2-$3\t/;  $h{$k}="$f\t$d1\t$d2\n"; }  END { if($pk && $h{$pk}) { print $h{$pk}; } }' | sort -k3,3nr > w${WIGGLE}.refseq_gencode.contained_not_just_overlapping

#3) CONTAINING (annotated within new exons)
#get ones which contain annotated
cat ${SPLIT}.raw | perl -ne 'BEGIN {$w='${WIGGLE}';} chomp; $f=$_; @f=split(/\t/,$f); ($s,$e)=($f[1],$f[2]); ($s2,$e2)=($f[8],$f[9]); $d1=abs($s2-$s); $d2=abs($e2-$e); $k="$c:$s-$e"; if($pk && $k ne $pk && $h{$pk}) { print $h{$pk}; delete $h{$pk}; } $pk=$k; if($s < $s2 && $e > $e2 && $d1 > $w && $d2 > $w) { $f=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/$1:$2-$3\t/;  $h{$k}="$f\t$d1\t$d2\n"; } END { if($pk && $h{$pk}) { print $h{$pk}; } }' | sort -k3,3nr > w${WIGGLE}.refseq_gencode.containing_annotated

#$BT window -a <(sort -k1,1 -k2,2n -k3,3n sim10.exons.split_ends.bed.w${WIGGLE}.raw.either_end) -b $ANNOT_EXONS -w $WIGGLE -v | uniq | perl -ne 'chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads)=split(/\t/,$f); print "$c:$s-$e\t$o\t$nr\t$reads\n";' >w${WIGGLE}.refseq_gencode.no_overlap_at_all
