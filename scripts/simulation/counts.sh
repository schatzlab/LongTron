#!/bin/bash
set -o pipefail -o nounset -o errexit
#count different slices of the output of comparison between 2 sets of exon coordinates (strand ignored)

simulated=$1
novel=$2
matching=$3
nmatching=$4

overlapping=$5
contained=$6
containing=$7
wiggle=$8

cut -f 1,2,3 $simulated | sort -u > ${simulated}.u
wc -l ${simulated}.u > ${simulated}.u.wc
total=`cut -d' ' -f 1 ${simulated}.u.wc`

cat /dev/null > all3
for f in $novel $matching $nmatching; do
	cut -f 1,2,3 $f | sort -u > ${f}.u
	wc -l ${f}.u > ${f}.u.wc
	cat ${f}.u >> all3
done
sort -u all3 > all3.u
all3=`wc -l all3.u | cut -d" " -f1`

perl -e '$w='${wiggle}'; $t='${total}'; $a='${all3}'; $n='`cut -d" " -f1 ${novel}.u.wc`'; chomp($n); $m='`cut -d" " -f1 ${matching}.u.wc`'; chomp($m); $nm='`cut -d" " -f1 ${nmatching}.u.wc`'; chomp($nm); if($t != $a) { print "total ".$t." DOES NOT match novel+matching+non-matching ".$a.", unexpected, quitting!\n"; } else { print "100%\t($t)\ttotal matches\n";} printf("breakdown of totals:\n%.1f%\t(%u)\tcompletely novel\n%.1f%\t(%u)\toverlap with both ends within wiggle ($w)\n%.1f%\t(%u)\toverlap but one or both ends not within wiggle ($w)\n",100*($n/$t),$n,100*($m/$t),$m,100*($nm/$t),$nm);'

cat /dev/null > all3.nmatching
for f in $overlapping $contained $containing; do
	cut -f 1 $f | sort -u | tr \\$ \\t > ${f}.u
	wc -l ${f}.u | cut -d' ' -f 1 > ${f}.u.wc
	cat ${f}.u >> all3.nmatching
done
sort -u all3.nmatching > all3.nmatching.u
all3=`wc -l all3.nmatching.u | cut -d" " -f1`

overlapping1=`wc -l ${overlapping}.u | cut -d" " -f1`
fgrep -v -f ${contained}.u ${overlapping}.u > ${overlapping}.u.not_contained
overlapping2=`wc -l ${overlapping}.u.not_contained | cut -d" " -f1`

containing1=`wc -l ${containing}.u | cut -d" " -f1`
fgrep -v -f ${contained}.u ${containing}.u > ${containing}.u.not_contained
containing2=`wc -l ${containing}.u.not_contained | cut -d" " -f1`

contained=`wc -l ${contained}.u | cut -d" " -f1`

perl -e '$w='${wiggle}'; $t='`cut -d" " -f1 ${nmatching}.u.wc`'; $a='${all3}'; if($t != $a) { print "total mismatches ".$t." DOES NOT match overlaps+contained+containing ".$a.", unexpected, quitting!\n"; } else { print "total non-matches ($t)\n";} $o='${overlapping2}'; chomp($o); $contained='${contained}'; chomp($contained); $containing='${containing2}'; chomp($containing); printf("breakdown of total non-matching:\n%.1f%\t(%u)\toverlapping annotated exons\n%.1f%\t(%u)\tcontaining annotated exons\n%.1f%\t(%u)\tcontained within annotated exons\n",100*($o/$t),$o,100*($containing/$t),$containing,100*($contained/$t),$contained);'

perl -e '$w='${wiggle}'; $t='${contained}'; print "total contained non-matches ($t)\n"; $o1='${overlapping1}'; chomp($o1); $o2='${overlapping2}'; chomp($o2); $containing1='${containing1}'; chomp($containing1); $containing2='${containing2}'; chomp($containing2); $contained_overlaps=$o1-$o2; $contained_containing=$containing1-$containing2; printf("breakdown of contained non-matches:\n%.1f%\t(%u)\toverlapping annotated exons\n%.1f%\t(%u)\tcontaining annotated exons\n",100*($contained_overlaps/$t),$contained_overlaps,100*($contained_containing/$t),$contained_containing);'

