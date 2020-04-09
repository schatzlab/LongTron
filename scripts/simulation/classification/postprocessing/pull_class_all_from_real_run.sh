#!/usr/bin/env bash
fl_training_examples_file=$1
#nonfl_training_examples_file=$1
nonfl_training_examples_file=$2

#this function keeps all the predictions where the predicted class (from the probabilities outut from the RF "5 class")
#matches the trained class
function find_matches_all  {
    f=$1
    start_pos=$2
    t=$3
    pos=`perl -e '%h=("problem-free"=>0,"non-recurrent"=>1,"recurrent"=>2,"novel"=>3); print "".($h{"'${f}'"}+'${start_pos}');'`
    echo "start_pos $start_pos f $f pos $pos"
    #filter first by whether the training category "class" (last column) matches what the current class we're compiling for (bash $f)
    #then if the training class (last column) matches the max prob. predicted class, then we assign a match ($m tracks this)
    #$c tracks the total number of training class assignments with predictions, whether those predictions are correct or not
    zcat all.${t}.overlaps.bgz | perl -ne 'BEGIN { $start='${start_pos}'-1; $pos='${pos}'-1; } chomp; $f=$_; @f=split(/\t/,$f); $class=pop(@f); next if($class ne "'${f}'"); $c++; $max=0; $max_k=-1; for($i=$start;$i<($start+4);$i++) { if($f[$i] > $max) { $max=$f[$i]; $max_k=$i; } } if($max_k == $pos) { print "$f\n"; $m++;} END { $a=$m/$c; print STDERR "'${f}'\t'${t}'\ttotal:$c\tmatches:$m\tratio:$a\n"; }' 2> ${f}.${t}.all.matches.info | sort -k${pos},${pos}nr > ${f}.${t}.all.matches
}

#find overlaps between 1) the predictions on the real data with 2) the training data with it's labels
#we want to see if predicted label and the training label are the same (keep the ones which are)
bedtools intersect -sorted -wo -a <(zcat features.classes.bgz) -b ${fl_training_examples_file}.all.bed.sorted | cut -f 1-3,7,44- | cut -f 1-16,20 | uniq | perl -ne 'chomp; $s=$_; $s=~s/\t/:/; $s=~s/\t/-/; print "$s\n";' | sort -u | bgzip > all.fl.overlaps.bgz
bedtools intersect -sorted -wo -a <(zcat features.classes.bgz) -b ${nonfl_training_examples_file}.all.bed.sorted | cut -f 1-3,7,44- | cut -f 1-16,20 | uniq | perl -ne 'chomp; $s=$_; $s=~s/\t/:/; $s=~s/\t/-/; print "$s\n";' | sort -u | bgzip > all.nonfl.overlaps.bgz

for f in problem-free non-recurrent recurrent novel; do
    #fl
    start_pos=3
    t="fl"
    find_matches_all $f $start_pos ${t}
    
    #nonfl
    start_pos=9
    t="nonfl"
    find_matches_all $f $start_pos ${t}
done

#then run to find, for each class (e.g. novel) the regions which overlap between fl and nonfl
#/bin/bash -x simulation/classification/postprocessing/find_fl_and_nonfl_agreement_matches.sh all.matches fl > find_fl_and_nonfl_agreement_matches.sh.run 2>&1

