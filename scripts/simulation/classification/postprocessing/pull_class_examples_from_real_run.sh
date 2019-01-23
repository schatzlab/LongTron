#!/usr/bin/env bash
fl_training_examples_file=$1
nonfl_training_examples_file=$1

function find_matches  {
    f=$1
    start_pos=$2
    t=$3
    n=$4
    pos=`perl -e '%h=("problem-free"=>0,"non-recurrent"=>1,"recurrent"=>2,"novel"=>3); print "".($h{"'${f}'"}+'${start_pos}');'`
    echo "start_pos $start_pos f $f pos $pos"
    cat ${f}.${t}.${n} | perl -ne 'BEGIN { $start='${start_pos}'-1; $pos='${pos}'-1; } chomp; $c++; $f=$_; @f=split(/\t/,$f); $max=0; $max_k=-1; for($i=$start;$i<($start+4);$i++) { if($f[$i] > $max) { $max=$f[$i]; $max_k=$i; } } if($max_k == $pos) { print "$f\n"; $m++;} END { $a=$m/$c; print STDERR "'${f}'\t'${t}'\ttotal:$c\tmatches:$m\tratio:$a\n"; }' > ${f}.${t}.${n}.matches 2> ${f}.${t}.${n}.matches.info
}

for f in problem-free non-recurrent recurrent novel; do
    #fl
    start_pos=3
    t="fl"
    coord=`egrep -e $f ${fl_training_examples_file}.head | cut -f 1`
    tabix features.classes.bgz $coord | cut -f 1-3,7,44- | perl -ne 'chomp; $s=$_; $s=~s/\t/:/; $s=~s/\t/-/; print "$s\n";' > ${f}.${t}.1
    find_matches $f $start_pos ${t} 1
    
    #nonfl
    start_pos=9
    t="nonfl"
    coord=`egrep -e $f ${nonfl_training_examples_file}.head | cut -f 1`
    tabix features.classes.bgz $coord | cut -f 1-3,7,44- | perl -ne 'chomp; $s=$_; $s=~s/\t/:/; $s=~s/\t/-/; print "$s\n";' > ${f}.${t}.1
    find_matches $f $start_pos ${t} 1
done

