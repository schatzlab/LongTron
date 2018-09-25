#!/bin/bash
set -o pipefail -o nounset -o errexit 

#base-0, so 5 runs if this is 4
#the following assumes that there are at least 3 runs
MAX_RUN_IDX=4
WIGGLE=20
ANNOTATION=gencode.v28.basic.annotation.transcripts_exon_count
JX_PREFIX='sim10.junctions.clean'

#get recurrent transcripts per categoy across all runs
echo "" > w${WIGGLE}.refseq_genecode.all_categories.trans.names.all
for f in `ls fl.${WIGGLE}.all.0/*.trans.names | cut -d'/' -f 2`; do
	comm -1 -2 fl.${WIGGLE}.all.0/$f fl.${WIGGLE}.all.1/$f > ${f}.all
	for i in $(seq `expr ${MAX_RUN_IDX} - 2` ${MAX_RUN_IDX}); do
		comm -1 -2 ${f}.all fl.${WIGGLE}.all.${i}/$f > ${f}.all.1
		mv ${f}.all.1 ${f}.all
	done
	cat ${f}.all >> w${WIGGLE}.refseq_genecode.all_categories.trans.names.all
	e=`wc -l ${f}.all | cut -d" " -f 1`
	#get total for this category across all runs
	wc -l fl.${WIGGLE}.all.?/${f} | egrep -v -e 'total' | perl -ne 'chomp; ($j,$n,$f)=split(/\s+/,$_); $s+=$n; $c++; END { $a=$s/$c; printf("'${f}': %.1f% (%u/%u)\n",100*('${e}'/$a),'${e}',$a); }'
	f2=`perl -e '$f1="'${f}'"; $f1=~s/\.names/\.counts/; print "$f1\n";'`
	fgrep -f ${f}.all fl.${WIGGLE}.all.?/$f2 | perl -ne 'chomp; ($f,$c,$n)=split(/\s+/,$_); $h{$n}+=$c; END { for $n (keys %h) { print "$n\t".$h{$n}."\n"; }}' | sort -k2,2nr > ${f}.all.counts
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

cat all.3.trans.counts | perl -ne 'BEGIN { open(IN,"<'${ANNOTATION}'"); %h; while($line=<IN>) { chomp($line); ($t,$c)=split(/\s+/,$line); $h{$t}=$c; } close(IN); } chomp; $f=$_; ($t,$c2)=split(/\t/,$f); $t=~s/\.[\+-]//; $c1=$h{$t}; $c3=$c2/$c1; print "$t\t$c3\n";' | sort -k2,2nr > all.3.trans.counts.normalized_by_exon_count

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
	cut -f 7 fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel | tr \; \\n | sort | uniq -c | sort -k1,1nr > fl.${WIGGLE}.all.${i}/${JX_PREFIX}.w${WIGGLE}.novel.names_counts
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
