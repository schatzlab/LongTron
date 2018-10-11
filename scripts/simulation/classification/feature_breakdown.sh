#!/bin/bash

#NOVELS=all.novels.names
NOVELS=recurrent.novels.names.sorted.u

echo -n 'category	' > all.tfeatures
cat ../transcript_stats/header >> all.tfeatures

#recurrent problems
~/filter_one_file_by_another.py -f all.3.trans.names.nonovels -t ../tfeatures -w -c 3 > all.3.trans.names.nonovels.tfeatures
cat all.3.trans.names.nonovels.tfeatures | perl -ne 'print "recurrent\t$_";' >> all.tfeatures

#non-problem, non-novels (good ones)
~/filter_one_file_by_another.py -f good_tnames.nonovels -t ../tfeatures -w -c 3 > good_tnames.nonovels.tfeatures
cat good_tnames.nonovels.tfeatures | perl -ne 'print "problem-free\t$_";' >> all.tfeatures

#novels
~/filter_one_file_by_another.py -f ${NOVELS} -t ../tfeatures -w -c 3 > ${NOVELS}.tfeatures
cat ${NOVELS}.tfeatures | perl -ne 'print "novel\t$_";' >> all.tfeatures

#non-recurrent problems
~/filter_one_file_by_another.py -f all.problem.names.sorted.u.nonall3.nonovels -t ../tfeatures -w -c 3 > all.problem.names.sorted.u.nonall3.nonovels.tfeatures
cat all.problem.names.sorted.u.nonall3.nonovels.tfeatures | perl -ne 'print "non-recurrent\t$_";' >> all.tfeatures

cat ../header.index | perl -ne 'chomp; ($i,$n)=split(/\t/,$_); $n=~s/\s+/_/g; $i++; `cut -f 1,$i all.tfeatures > all.tfeatures.$n`;'
