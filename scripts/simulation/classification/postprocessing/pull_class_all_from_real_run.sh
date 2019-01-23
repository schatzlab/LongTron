#!/usr/bin/env bash
fl_training_examples_file=$1
nonfl_training_examples_file=$1

function find_matches_all  {
    f=$1
    start_pos=$2
    t=$3
    pos=`perl -e '%h=("problem-free"=>0,"non-recurrent"=>1,"recurrent"=>2,"novel"=>3); print "".($h{"'${f}'"}+'${start_pos}');'`
    echo "start_pos $start_pos f $f pos $pos"
    zcat all.${t}.overlaps.bgz | perl -ne 'BEGIN { $start='${start_pos}'-1; $pos='${pos}'-1; } chomp; $f=$_; @f=split(/\t/,$f); $class=pop(@f); next if($class ne "'${f}'"); $c++; $max=0; $max_k=-1; for($i=$start;$i<($start+4);$i++) { if($f[$i] > $max) { $max=$f[$i]; $max_k=$i; } } if($max_k == $pos) { print "$f\n"; $m++;} END { $a=$m/$c; print STDERR "'${f}'\t'${t}'\ttotal:$c\tmatches:$m\tratio:$a\n"; }' > ${f}.${t}.all.matches 2> ${f}.${t}.all.matches.info
}

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
