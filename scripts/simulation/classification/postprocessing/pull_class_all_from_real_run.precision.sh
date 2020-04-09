#!/usr/bin/env bash
#find the precision for each class for both FL and nonFL types
#THIS ASSUMES  pull_class_all_from_real_run.sh was run first to get the 1) overlap files and 2) to get the number of true predictions
#both of which are used below

recall_table=$1

#this function keeps all the predictions where the predicted class (from the probabilities outut from the RF "5 class")
#matches the trained class
function find_matches_all  {
    f=$1
    start_pos=$2
    t=$3
    pos=`perl -e '%h=("problem-free"=>0,"non-recurrent"=>1,"recurrent"=>2,"novel"=>3); print "".($h{"'${f}'"}+'${start_pos}');'`
    echo "start_pos $start_pos f $f pos $pos"
    #want to find the total number of predictions for the current class ($pred_class)
    #we'll get the correct ones from the recall numbers calculated separately
    zcat all.${t}.overlaps.bgz | perl -ne 'BEGIN { $pred_class="'$f'"; $start='${start_pos}'-1; $pos='${pos}'-1; } chomp; $f=$_; @f=split(/\t/,$f); $class=pop(@f); $c++; $max=0; $max_k=-1; for($i=$start;$i<($start+4);$i++) { if($f[$i] > $max) { $max=$f[$i]; $max_k=$i; } } if($max_k == $pos) { print "$f\n"; $m++;} END { $a=$m/$c; print STDERR "'${f}'\t'${t}'\t$m\t$c\n"; }' 2>> all.predictions.info | sort -k${pos},${pos}nr > ${f}.${t}.all.predictions
}

echo -n "" > all.predictions.info
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

#now join # of predictions with # of true predictions to get precision
cat <(cut -f 1-3 all.predictions.info) $recall_table | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); $k=$f[0]."\t".$f[1]; if(scalar(@f) == 3) { $h{$k}=$f[2]; next; } $total_preds = $h{$k}; $tp = $f[3]; $tp=~s/matches://; $r=$tp/$total_preds; print "$k\ttotal:$total_preds\tmatches:$tp\tratio:$r\n";' > all.recall.info
