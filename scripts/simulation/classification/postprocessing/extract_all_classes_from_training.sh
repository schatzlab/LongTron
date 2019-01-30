#!/usr/bin/env bash

#training features file (used as input to the training)
features_f=$1
#fl or nonfl
training_type=$2

#the following is just for a few examples
#echo -n "" > all_classes.head
#for f in problem-free non-recurrent recurrent novel; do
#    egrep -e '	'$f'	' $features_f | head | perl -ne 'chomp; $s=$_; $s=~s/\t/:/; $s=~s/\t/-/; print "$s\n";' > ${f}.head
#    head -1 ${f}.head | cut -f 1 | perl -ne 'chomp; $f=$_; print "$f\t'${f}'\n";' >> all_classes.head
#done

#this does all the data
echo -n "" > ${run_type}.all_classes.all.bed
for f in problem-free non-recurrent recurrent novel; do
    #egrep -e '	'$f'	' $features_f | perl -ne 'chomp; $s=$_; $s=~s/\t/:/; $s=~s/\t/-/; print "$s\n";' | cut -f1,6 >> all_classes.all
    #training data is already 0-based coord so fits BED format
    egrep -e '	'$f'	' $features_f | cut -f1-3,8 >> ${run_type}.all_classes.all.bed
done
sort -k1,1 -k2,2n -k3,3n ${run_type}.all_classes.all.bed > ${run_type}.all_classes.all.bed.sorted
rm ${run_type}.all_classes.all.bed
