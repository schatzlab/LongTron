#!/bin/bash
#generate a set of 1) training 2) validation and 3) test partitions for each category of transcript/reads

ntraining=$1
nvalidation=$2
input=$3

cut -f 8- $input > ${input}.x
#cut -f 1 ${input}.x | perl -ne 'BEGIN { %a=("problem-free"=>1,"non-recurrent"=>2,"recurrent"=>3,"novel"=>4); } chomp; $r=$_; @b=(0,0,0,0); $b[$a{$r}-1]=1; print "$r\t".join("\t",@b)."\n";' > ${input}.y
#redo 4-class with 0,1,2,3 class labels
cut -f 1 ${input}.x | perl -ne 'BEGIN { %a=("problem-free"=>0,"non-recurrent"=>1,"recurrent"=>2,"novel"=>3); } chomp; $r=$_; print "$r\t".$a{$r}."\n";' > ${input}.y
ln -fs ${input}.x features.x
ln -fs ${input}.y features.y

echo -n "" > training.x
echo -n "" > training.y
echo -n "" > validation.x
echo -n "" > validation.y
echo -n "" > testing.x
echo -n "" > testing.y
echo -n "" > num_training_samples_per_cat
for f in "problem-free" "non-recurrent" "recurrent" "novel";
do
	#first we grep out the specific set we want and shuffle (randomize) the lines
	egrep -e "^${f}" features.x | cut -f2- | shuf > features.x.${f}
	egrep -e "^${f}" features.y | cut -f2- > features.y.${f}
	wc -l features.y.${f} >> num_training_samples_per_cat

	#now figure out sizes based on fractions passed in for 1) training and 2) validation sets
	size=`wc -l features.x.${f} | cut -d" " -f1`

	#training first
	t=`perl -e '$total='${size}'; $frac='${ntraining}'; $size=int($frac*$total); print "$size\n";'`
	head -${t} features.x.${f} >> training.x
	head -${t} features.y.${f} >> training.y

	#validation middle
	v=`perl -e '$total='${size}'; $frac='${nvalidation}'; $size=int($frac*$total); print "$size\n";'`
	#need to add +1 to the original training size so we can start on the next line
	t2=`perl -e '$total='${size}'; $frac='${ntraining}'; $size=int($frac*$total)+1; print "$size\n";'`
	tail -n+${t2} features.x.${f} | head -${v} >> validation.x	
	head -${v} features.y.${f} >> validation.y

	#finally testing
	vt=`perl -e '$total='${size}'; $t='${t}'; $v='${v}'; $size=$total-($t+$v); print "$size\n";'`
	tail -n ${vt} features.x.${f} >> testing.x
	head -${v} features.y.${f} >> testing.y
done

#now do a balanced version of the training set
echo -n "" > training.x.bal
echo -n "" > training.y.bal
min=`cat num_training_samples_per_cat | perl -ne 'BEGIN { $min=(2**32)-1; } chomp; ($c,$name)=split(/\s+/,$_); next if($c > $min); $min=$c; END { print "$min\n";}'`
for f in "problem-free" "non-recurrent" "recurrent" "novel";
do
	head -${min} features.x.${f} >> training.x.bal
	head -${min} features.y.${f} >> training.y.bal
done


#create 1) binary (0,1) labels and separately 2) only use recurrent and problem-free
n=`head -1 training.x | tr \\\t \\\n | wc -l | cut -d" " -f 1`
for t in "training" "validation" "testing";
do
	#cat ${t}.y | perl -ne 'chomp; ($problem_free,$non_recurrent,$recurrent,$novel) = split(/\t/,$_); if($problem_free == 1) { print "0\n"; next; } print "1\n";' > binary.y.${t}
	cat ${t}.y | perl -ne 'chomp; $f=$_; if($f != 0) { print "1\n"; next; } print "0\n";' > binary.y.${t}
	#binary: recurrent vs. problem-free
	#paste ${t}.x ${t}.y | egrep -e '	.	0	.	0$' | perl -ne 'chomp; @f=split(/\t/,$_); @f1=splice(@f,0,'${n}'); print "".join("\t",@f1)."\n"; ($problem_free,$non_recurrent,$recurrent,$novel) = @f; print STDERR "$recurrent\n";' > ${t}.2.x 2> ${t}.2.y
	#binary: recurrent + non-recurrent vs. problem-free
	#paste ${t}.x ${t}.y | egrep -e '	.	.	.	0$' | perl -ne 'chomp; @f=split(/\t/,$_); @f1=splice(@f,0,'${n}'); print "".join("\t",@f1)."\n"; ($problem_free,$non_recurrent,$recurrent,$novel) = @f; print STDERR "$problem_free\n";' > ${t}.3a.x 2> ${t}.3a.y
	#binary: recurrent + novels vs. problem-free
	#paste ${t}.x ${t}.y | egrep -e '	.	0	.	.$' | perl -ne 'chomp; @f=split(/\t/,$_); @f1=splice(@f,0,'${n}'); print "".join("\t",@f1)."\n"; ($problem_free,$non_recurrent,$recurrent,$novel) = @f; print STDERR "$problem_free\n";' > ${t}.3b.x 2> ${t}.3b.y
done

