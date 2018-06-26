#!/bin/bash

echo "reads_required	novel_junctions	in_pacbio" > original_plus_4_wiggle.juncs.sorted.full.novel
for i in 1 2 5 7 10 15 20 50 100 500 1000
do
	echo -n "$i	" >> original_plus_4_wiggle.juncs.sorted.full.novel 
	#cat original_plus_5_wiggle.juncs.sorted.full | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); if($f[4] >= '${i}' && $f[9]==0 && $f[14]==0 && $f[19]==0 && $f[24]==0 && $f[29]==0) { print "$f\n";}' | cut -f1-3 > original_plus_5_wiggle.juncs.sorted.full.novel_${i} 
	cat original_plus_5_wiggle.juncs.sorted.full | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); if($f[4] >= '${i}' && $f[9]==0 && $f[14]==0 && $f[19]==0 && $f[24]==0) { $p=0; $p=1 if($f[29] > 0); print "".join("\t",($f[0],$f[1],$f[2]))."\t$p\n"; }' > original_plus_4_wiggle.juncs.sorted.full.novel_${i} 
	cat original_plus_4_wiggle.juncs.sorted.full.novel_${i} | wc -l >> original_plus_4_wiggle.juncs.sorted.full.novel
done
