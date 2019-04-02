#!/bin/bash
set -o pipefail -o errexit 

FILTER=~/filter_one_file_by_another.py

#e.g. gencode.v28.basic.annotation.transcripts_exon_count
ANNOTATION=$1
FIX_TRANSCRIPT_NAMES=0
if [ -n "$2" ]; then
	FIX_TRANSCRIPT_NAMES=$2
fi
#base-0, so 5 runs if this is 4
#the following assumes that there are at least 3 runs
MAX_RUN_IDX=4
WIGGLE=20
JX_PREFIX='sim10.junctions.clean'

#get recurrent transcripts per categoy across all runs
echo -n "" > w${WIGGLE}.refseq_genecode.all_categories.trans.names.all
for f in `ls fl.${WIGGLE}.all.0/*.trans.names | cut -d'/' -f 2`; do
	#need to trim transcript names first
	if [ ! -f fl.${WIGGLE}.all.0/${f}.pre_trimming ]; then
		cp fl.${WIGGLE}.all.0/${f} fl.${WIGGLE}.all.0/${f}.pre_trimming
	fi
	cat fl.${WIGGLE}.all.0/${f}.pre_trimming | perl -ne 'chomp; $t=$_; if('${FIX_TRANSCRIPT_NAMES}'==1 && $t !~ /PAR/) { $t=~s/^([^\.]+\.\d+).*$/$1/; } print "$t\n";' > fl.${WIGGLE}.all.0/${f}
	if [ ! -f fl.${WIGGLE}.all.1/${f}.pre_trimming ]; then
		cp fl.${WIGGLE}.all.1/${f} fl.${WIGGLE}.all.1/${f}.pre_trimming
	fi
	cat fl.${WIGGLE}.all.1/${f}.pre_trimming | perl -ne 'chomp; $t=$_; if('${FIX_TRANSCRIPT_NAMES}'==1 && $t !~ /PAR/) { $t=~s/^([^\.]+\.\d+).*$/$1/; } print "$t\n";' > fl.${WIGGLE}.all.1/${f}
	#now start by comparing
	comm -1 -2 fl.${WIGGLE}.all.0/$f fl.${WIGGLE}.all.1/$f > ${f}.all
	for i in $(seq `expr ${MAX_RUN_IDX} - 2` ${MAX_RUN_IDX}); do
		if [ ! -f fl.${WIGGLE}.all.${i}/${f}.pre_trimming ]; then
			cp fl.${WIGGLE}.all.${i}/${f} fl.${WIGGLE}.all.${i}/${f}.pre_trimming
		fi
		cat fl.${WIGGLE}.all.${i}/${f}.pre_trimming | perl -ne 'chomp; $t=$_; if('${FIX_TRANSCRIPT_NAMES}'==1 && $t !~ /PAR/) { $t=~s/^([^\.]+\.\d+).*$/$1/; } print "$t\n";' > fl.${WIGGLE}.all.${i}/${f}
		comm -1 -2 ${f}.all fl.${WIGGLE}.all.${i}/$f > ${f}.all.1
		mv ${f}.all.1 ${f}.all
	done
	cat ${f}.all >> w${WIGGLE}.refseq_genecode.all_categories.trans.names.all
	e=`wc -l ${f}.all | cut -d" " -f 1`
	#get total for this category across all runs
	wc -l fl.${WIGGLE}.all.?/${f} | egrep -v -e 'total' | perl -ne 'chomp; ($j,$n,$f)=split(/\s+/,$_); $s+=$n; $c++; END { $a=$s/$c; printf("'${f}': %.1f% (%u/%u)\n",100*('${e}'/$a),'${e}',$a); }'
	f2=`perl -e '$f1="'${f}'"; $f1=~s/\.names/\.counts/; print "$f1\n";'`
	#fgrep -f ${f}.all fl.${WIGGLE}.all.?/$f2 | perl -ne 'chomp; ($f,$c,$n)=split(/\s+/,$_); if('${FIX_TRANSCRIPT_NAMES}'==1) { $n=~s/^([^\.]+\.\d+).*$/$1/; } $h{$n}+=$c; END { for $n (keys %h) { print "$n\t".$h{$n}."\n"; }}' | sort -k2,2nr > ${f}.all.counts
	$FILTER -f ${f}.all -t fl.${WIGGLE}.all.?/$f2 -w -c 2 -p'\s+' | perl -ne 'chomp; ($f,$c,$n)=split(/\s+/,$_); if('${FIX_TRANSCRIPT_NAMES}'==1 && $n !~ /PAR/) { $n=~s/^([^\.]+\.\d+).*$/$1/; } $h{$n}+=$c; END { for $n (keys %h) { print "$n\t".$h{$n}."\n"; }}' | sort -k2,2nr > ${f}.all.counts
done
sort -u w${WIGGLE}.refseq_genecode.all_categories.trans.names.all > w${WIGGLE}.refseq_genecode.all_categories.trans.names.all.sorted

#now get ones shared across all 3 categories across all runs
if [ -f all.3.trans.names ];
then
	rm all.3.trans.names
fi
for f in `ls fl.${WIGGLE}.all.0/*.trans.names | cut -d'/' -f 2`; do
	if [[ ! -f all.3.trans.names ]]; then
		cat ${f}.all | sort > all.3.trans.names
	fi
	comm -1 -2 all.3.trans.names <(sort ${f}.all) | sort > all.3.trans.names.1
	mv all.3.trans.names.1 all.3.trans.names
done
fgrep -f all.3.trans.names *.all.counts | perl -ne 'chomp; ($j,$n)=split(/:/,$_); ($n,$c)=split(/\s+/,$n); $h{$n}+=$c; END { for $n (keys %h) { print "$n\t".$h{$n}."\n"; }}' | sort -k2,2nr > all.3.trans.counts

cat all.3.trans.counts | perl -ne 'BEGIN { open(IN,"<'${ANNOTATION}'"); %h; while($line=<IN>) { chomp($line); ($t,$c)=split(/\s+/,$line); $h{$t}=$c; } close(IN); } chomp; $f=$_; ($t,$c2)=split(/\t/,$f); if('${FIX_TRANSCRIPT_NAMES}'==1 && $t !~ /PAR/) { $t=~s/^([^\.]+\.\d+).*$/$1/; } $c1=$h{$t}; $c3=$c2/$c1; print "$t\t$c3\n";' | sort -k2,2nr > all.3.trans.counts.normalized_by_exon_count

echo "" > all.problem.names
for i in $(seq 0 ${MAX_RUN_IDX}); do cat fl.${WIGGLE}.all.${i}/*.names | cut -d'_' -f 1 | cut -d'.' -f 1,2 >> all.problem.names ; done
sort -u all.problem.names | egrep -v -e '^$' > all.problem.names.sorted.u

#novel names
echo "" > all.novels.names
if [ -f recurrent.novels.names ]
then
	rm recurrent.novels.names
fi
for i in $(seq 0 ${MAX_RUN_IDX});
do
	cut -f 7 fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel | tr \; \\n | cut -d'_' -f 1 | cut -d'.' -f 1,2 | sort | uniq -c | sort -k1,1nr > fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel.names_counts
	cat fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel.names_counts | perl -ne 'chomp; ($j,$c,$n)=split(/\s+/,$_); print "$n\n";' | sort > fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel.names
	cat fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel.names >> all.novels.names
	if [ ! -f recurrent.novels.names ]
	then
		cat fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel.names > recurrent.novels.names
	fi
	comm -1 -2 fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel.names recurrent.novels.names | sort > recurrent.novels.names.1
	mv recurrent.novels.names.1 recurrent.novels.names
done
sort -u recurrent.novels.names | egrep -v -e '^$' > recurrent.novels.names.sorted.u

NOVELS=recurrent.novels.names.sorted.u
ALL_CLASSES=class2transcript.all

#novels
cat $NOVELS | perl -ne 'chomp; print "novel\t$_\n";' > $ALL_CLASSES
#recurrent
$FILTER -f ${NOVELS} -t all.3.trans.names -w -c 0 -n > all.3.trans.names.nonovels
cat all.3.trans.names.nonovels | perl -ne 'chomp; print "recurrent\t$_\n";' >> $ALL_CLASSES
#non-recurrent
$FILTER -f all.3.trans.names -t all.problem.names.sorted.u -w -c 0 -n > all.problem.names.sorted.u.nonall3
$FILTER -f ${NOVELS} -t all.problem.names.sorted.u.nonall3 -w -c 0 -n > all.problem.names.sorted.u.nonall3.nonovels
cat all.problem.names.sorted.u.nonall3.nonovels | perl -ne 'chomp; print "non-recurrent\t$_\n";' >> $ALL_CLASSES
#problem-free
$FILTER -f all.problem.names.sorted.u -t ${ANNOTATION} -w -c 0 -n | cut -f 1 > good_tnames
$FILTER -f ${NOVELS} -t good_tnames -w -c 0 -n > good_tnames.nonovels
cat good_tnames.nonovels | perl -ne 'chomp; print "problem-free\t$_\n";' >> $ALL_CLASSES
