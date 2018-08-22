#!/bin/bash
set -o pipefail -o nounset -o errexit 

for f in `ls fl.20.all.1/*.trans.names | cut -d'/' -f 2`; do
	comm -1 -2 fl.20.all.0/$f fl.20.all.1/$f > ${f}.all
	for i in 2 3 4; do
		comm -1 -2 ${f}.all fl.20.all.${i}/$f > ${f}.all.1
		mv ${f}.all.1 ${f}.all
	done
	e=`wc -l ${f}.all | cut -d" " -f 1`
	wc -l fl.20.all.?/${f} | egrep -v -e 'total' | perl -ne 'chomp; ($j,$n,$f)=split(/\s+/,$_); $s+=$n; $c++; END { $a=$s/$c; printf("'${f}': %.1f% (%u/%u)\n",100*('${e}'/$a),'${e}',$a); }'
	f2=`perl -e '$f1="'${f}'"; $f1=~s/\.names/\.counts/; print "$f1\n";'`
	fgrep -f ${f}.all fl.20.all.?/$f2 | perl -ne 'chomp; ($f,$c,$n)=split(/\s+/,$_); $h{$n}+=$c; END { for $n (keys %h) { print "$n\t".$h{$n}."\n"; }}' | sort -k2,2nr > ${f}.all.counts
done


for f in `ls fl.20.all.1/*.trans.names | cut -d'/' -f 2`; do
	if [[ ! -f all.3.trans.names ]]; then
		cat ${f}.all > all.3.trans.names
	fi
	comm -1 -2 all.3.trans.names ${f}.all > all.3.trans.names.1
	mv all.3.trans.names.1 all.3.trans.names
done
fgrep -f all.3.trans.names *.all.counts | perl -ne 'chomp; ($j,$n)=split(/:/,$_); ($n,$c)=split(/\s+/,$n); $h{$n}+=$c; END { for $n (keys %h) { print "$n\t".$h{$n}."\n"; }}' | sort -k2,2nr > all.3.trans.counts

