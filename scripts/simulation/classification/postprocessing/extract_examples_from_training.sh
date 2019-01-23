#!/usr/bin/env bash
features_f=$1

echo -n "" > all_classes.head
for f in problem-free non-recurrent recurrent novel; do
    egrep -e '	'$f'	' $features_f | head | perl -ne 'chomp; $s=$_; $s=~s/\t/:/; $s=~s/\t/-/; print "$s\n";' > ${f}.head
    head -1 ${f}.head | cut -f 1 | perl -ne 'chomp; $f=$_; print "$f\t'${f}'\n";' >> all_classes.head
done

echo -n "" > all_classes.all.bed
for f in problem-free non-recurrent recurrent novel; do
    #egrep -e '	'$f'	' $features_f | perl -ne 'chomp; $s=$_; $s=~s/\t/:/; $s=~s/\t/-/; print "$s\n";' | cut -f1,6 >> all_classes.all
    #training data is already 0-based coord so fits BED format
    egrep -e '	'$f'	' $features_f | cut -f1-3,8 >> all_classes.all.bed
done
sort -k1,1 -k2,2n -k3,3n all_classes.all.bed > all_classes.all.bed.sorted
rm all_classes.all.bed
