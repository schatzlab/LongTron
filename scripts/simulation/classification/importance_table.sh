#!/usr/bin/env bash

echo 'Category	First	Second	Third	Fourth	Fifth'

for f in Oxford_FL.txt Oxford_nFL.txt PacBio_FL.txt PacBio_nFL.txt ; do
    name=`perl -e '$n="'$f'"; $n=~s/\.txt$//; print "$n";'`
    grep importance $f | head -1 | tr -s ' '  \\n | fgrep -v importance | tr : \\t | perl -pe 'print "". ++$i.":";' | join -t'	' -j 1 - feature_map.txt | sort -t'	' -k2,2gr > ${f}_4_importances.all.tsv
    head -5 ${f}_4_importances.all.tsv | perl -ne 'BEGIN { $n="'$name'"; $n=~s/_/ /; print "$n 4"; } chomp; @f=split(/\t/,$_); printf("\t%s (%.0f\%)",$f[2],100*$f[1]); END { print "\n";}'
    grep importance $f | tail -n1 | tr -s ' '  \\n | fgrep -v importance | tr : \\t | perl -pe 'print "". ++$i.":";' | join -t'	' -j 1 - feature_map.txt | sort -t'	' -k2,2gr > ${f}_2_importances.all.tsv
    head -5 ${f}_2_importances.all.tsv | perl -ne 'BEGIN { $n="'$name'"; $n=~s/_/ /; print "$n 2"; } chomp; @f=split(/\t/,$_); printf("\t%s (%.0f\%)",$f[2],100*$f[1]); END { print "\n";}'
done
