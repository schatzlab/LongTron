#!/bin/bash
set -o pipefail -o nounset -o errexit 

scripts=`perl -e '$f="'${0}'"; $f=~s/\/[^\/]+$/\//; print "$f\n";'`

BT=/data/bedtools2/bin/bedtools

JUNCTIONS=$1
#refseq_junctions1.nostrand.split_ends.bed
#gencode_refseq_junctions.nostrand.sorted
#gencode.v28.basic.annotation.junctions.sorted
ANNOT_JUNCTIONS=$2
WIGGLE=$3
TID_FILE=$4

#filter out any non-cannoical mapped junction
egrep -v -e '(HLA)|(alt)|(random)|(decoy)|(chrUn)' ${JUNCTIONS} > ${JUNCTIONS}.clean

#BED format
cat ${JUNCTIONS}.clean | perl -ne 'chomp; ($c,$s,$e,$o,$nr,$reads)=split(/\t/,$_); $s--; print "".join("\t",($c,$s,$e,$o,$nr,$reads))."\n";' > ${JUNCTIONS}.clean.bed
cat $ANNOT_JUNCTIONS | perl -ne 'chomp; ($c,$s,$e,$t,$eid)=split(/\t/,$_); $s--; print "$c\t$s\t$e\t$t\t$eid\n";' > ${ANNOT_JUNCTIONS}.bed


#totally novel with wiggle (no split ends complications)
$BT window -a ${JUNCTIONS}.clean.bed -b ${ANNOT_JUNCTIONS}.bed -w $WIGGLE -v | perl -ne 'BEGIN { $w='${WIGGLE}'; open(IN,"<'${TID_FILE}'"); %tids; while($line=<IN>) { chomp($line); ($tid,$n)=split(/\t/,$line); $tids{$tid}=$n; } close(IN); } chomp; $f=$_; ($c1,$s1,$e1,$o,$nr,$reads)=split(/\t/,$f); %mtids; for $tid (split(/,/,$reads)) { ($tid,$eid)=split(/\./,$tid); $mtids{$tids{$tid}}=1; } print "$c1\t$s1\t$e1\t$o\t$nr\t$reads\t".(join(";",sort keys %mtids))."\n";' > sim10.junctions.clean.w${WIGGLE}.novel

#determine matching/nonmatching of junctions assigned to reads/transcripts
#To match:
#1) each end must be within window of it's annotated end, this decides which file it goes into (doesn't depend on #2)
#2) at least one  of the aligned jx's reads must match the annotated transcript ID it's overlapping
#and the exon/jx idx must match within the overlapping transcript, this decides if it gets a 1 or 0 for the matching column in output


$BT window -a ${JUNCTIONS}.clean.bed -b ${ANNOT_JUNCTIONS}.bed -w $WIGGLE | perl -ne 'BEGIN { $w='${WIGGLE}'; open(IN,"<'${TID_FILE}'"); %tids; while($line=<IN>) { chomp($line); ($tid,$n)=split(/\t/,$line); $n=~s/_0_[+-]$//; $tids{$tid}=uc($n); } close(IN); %mtids; } chomp; $f=$_; ($c1,$s1,$e1,$o,$nr,$reads,$c2,$s2,$e2,$t,$eid)=split(/\t/,$f); $t=~s/_0_[+-]$//; $t=uc($t); $d1=abs($s1-$s2); $d2=abs($e1-$e2); $c3=($d1>$d2?$s1:$e1); $k="$c1\t$s1\t$e1"; $matched_tid=0; for $tid (split(/,/,$reads)) { ($tid,$eid2)=split(/\./,$tid); if($tids{$tid} eq $t && $eid eq $eid2 && $mtids{$k} != 2) { $matched_tid=1; $mtids{$k}=1; break; } }  $line=join("\t",($c1,$s1,$e1,$o,$nr,$reads,$c3)); if($d1 <= $w && $d2 <= $w) { $mtids{$k} = 2 if($matched_tid == 1); delete $h2{$k}; $h{$k}=$line; } else { if(!$h{$k}) { $h2{$k}=$line; } }  END { for $k (keys %h2) { $matched=$mtids{$k}; $matched=0 if(!$mtids{$k}); $k=$h2{$k}; print "$k\t$matched\n"; } for $k (keys %h) { $matched=$mtids{$k}; $matched=0 if(!$mtids{$k}); $k=$h{$k}; print STDERR "$k\t$matched\n"; }}' 2> sim10.junctions.clean.w${WIGGLE}.matching > sim10.junctions.clean.w${WIGGLE}.nonmatching 

ALL_NONMATCH=sim10.junctions.split_ends.bed.w${WIGGLE}.either_end.sorted
cat sim10.junctions.clean.w${WIGGLE}.nonmatching | sort -k1,1 -k2,2n -k3,3n > $ALL_NONMATCH

#now separate into further categories those which have some form of overlap:

$BT window -a $ALL_NONMATCH -b ${ANNOT_JUNCTIONS}.bed -w $WIGGLE | uniq > ${ALL_NONMATCH}.raw

#1) OVERLAPPING
#use the original full coord range (not split ends) for both of the sets of intervals this time to remove contained/overlapping junctions
cat ${ALL_NONMATCH}.raw | perl -ne 'BEGIN {$w='${WIGGLE}';} chomp; $f=$_; @f=split(/\t/,$f); ($s,$e)=($f[1],$f[2]); ($s2,$e2)=($f[9],$f[10]); $d1=abs($s2-$s); $d2=abs($e2-$e); $k="$c:$s-$e"; if($pk && $k ne $pk && $h{$pk}) { print $h{$pk}; delete $h{$pk}; } $pk=$k; if( ($s < $s2 && $d1 > $w  && !($e > $e2 && $d2 > $w)) || ($e > $e2 && $d2 > $w && !($s < $s2 && $d1 > $w)) ) { $f=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/$1:$2-$3\t/;  $h{$k}="$f\t$d1\t$d2\n"; } END { if($pk && $h{$pk}) { print $h{$pk}; } }' | sort -k3,3nr > w${WIGGLE}.refseq_gencode.overlapping_not_contained

#2) CONTAINED (new junctions within annotated)
#get the contained ones
cat ${ALL_NONMATCH}.raw | perl -ne 'BEGIN {$w='${WIGGLE}';} chomp; $f=$_; @f=split(/\t/,$f); ($s,$e)=($f[1],$f[2]); ($s2,$e2)=($f[9],$f[10]); $d1=abs($s2-$s); $d2=abs($e2-$e); $k="$c:$s-$e"; if($pk && $k ne $pk && $h{$pk}) { print $h{$pk}; delete $h{$pk}; } $pk=$k; if( ($s >= $s2 && $e <= $e2) || ($s < $s2 && $d1 <= $w && !($e > $e2 && $d2 > $w)) || ($e > $e2 && $d2 <= $w && !($s < $s2 && $d1 > $w)) ) { $f=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/$1:$2-$3\t/;  $h{$k}="$f\t$d1\t$d2\n"; }  END { if($pk && $h{$pk}) { print $h{$pk}; } }' | sort -k3,3nr > w${WIGGLE}.refseq_gencode.contained_not_just_overlapping

#3) CONTAINING (annotated within new junctions)
#get ones which contain annotated
cat ${ALL_NONMATCH}.raw | perl -ne 'BEGIN {$w='${WIGGLE}';} chomp; $f=$_; @f=split(/\t/,$f); ($s,$e)=($f[1],$f[2]); ($s2,$e2)=($f[9],$f[10]); $d1=abs($s2-$s); $d2=abs($e2-$e); $k="$c:$s-$e"; if($pk && $k ne $pk && $h{$pk}) { print $h{$pk}; delete $h{$pk}; } $pk=$k; if($s < $s2 && $e > $e2 && $d1 > $w && $d2 > $w) { $f=~s/^([^\t]+)\t([^\t]+)\t([^\t]+)\t/$1:$2-$3\t/;  $h{$k}="$f\t$d1\t$d2\n"; } END { if($pk && $h{$pk}) { print $h{$pk}; } }' | sort -k3,3nr > w${WIGGLE}.refseq_gencode.containing_annotated

/bin/bash -x $scripts/counts.sh ${JUNCTIONS}.clean.bed sim10.junctions.clean.w${WIGGLE}.novel sim10.junctions.clean.w${WIGGLE}.matching sim10.junctions.clean.w${WIGGLE}.nonmatching w${WIGGLE}.refseq_gencode.overlapping_not_contained w${WIGGLE}.refseq_gencode.contained_not_just_overlapping w${WIGGLE}.refseq_gencode.containing_annotated ${WIGGLE} ${TID_FILE} > ${1}.w${WIGGLE}.counts.tsv

ln -fs ${1}.w${WIGGLE}.counts.tsv counts.out
