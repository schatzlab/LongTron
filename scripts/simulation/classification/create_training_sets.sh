#!/bin/bash
#generate a set of 1) training 2) validation and 3) test partitions for each category of transcript/reads

ntraining=$1
nvalidation=$2

cut -f 1,11,12,14,18,19,21- all.tfeatures | tail -n+2 | perl -ne 'BEGIN { %h=("non-recurrent"=>0,"novel"=>1,"problem-free"=>2,"recurrent"=>3); } chomp; @f=split(/\t/,$_); $n=shift(@f); @a=(0,0,0,0); $a[$h{$n}]=1; print "$n\t".join("\t",@f)."\n"; print STDERR "$n\t".join("\t",@a)."\n";' > all.tfeatures.all_x_vectors 2> all.tfeatures.all_y_vectors

echo -n "" > all.tfeatures.all_x_vectors.training
echo -n "" > all.tfeatures.all_y_vectors.training
echo -n "" > all.tfeatures.all_x_vectors.validation
echo -n "" > all.tfeatures.all_y_vectors.validation
echo -n "" > all.tfeatures.all_x_vectors.testing
echo -n "" > all.tfeatures.all_y_vectors.testing
for f in "non-recurrent" "novel" "problem-free" "recurrent";
do
	#first we grep out the specific set we want and shuffle (randomize) the lines
	egrep -e "^${f}" all.tfeatures.all_x_vectors | cut -f2- | shuf > all.tfeatures.all_x_vectors.${f}
	egrep -e "^${f}" all.tfeatures.all_y_vectors | cut -f2- > all.tfeatures.all_y_vectors.${f}

	#now figure out sizes based on fractions passed in for 1) training and 2) validation sets
	size=`wc -l all.tfeatures.all_x_vectors.${f} | cut -d" " -f1`

	#training first
	t=`perl -e '$total='${size}'; $frac='${ntraining}'; $size=int($frac*$total); print "$size\n";'`
	head -${t} all.tfeatures.all_x_vectors.${f} >> all.tfeatures.all_x_vectors.training
	head -${t} all.tfeatures.all_y_vectors.${f} >> all.tfeatures.all_y_vectors.training

	#validation middle
	v=`perl -e '$total='${size}'; $frac='${nvalidation}'; $size=int($frac*$total); print "$size\n";'`
	#need to add +1 to the original training size so we can start on the next line
	t2=`perl -e '$total='${size}'; $frac='${ntraining}'; $size=int($frac*$total)+1; print "$size\n";'`
	tail -n+${t2} all.tfeatures.all_x_vectors.${f} | head -${v} >> all.tfeatures.all_x_vectors.validation	
	head -${v} all.tfeatures.all_y_vectors.${f} >> all.tfeatures.all_y_vectors.validation

	#finally testing
	vt=`perl -e '$total='${size}'; $t='${t}'; $v='${v}'; $size=$total-($t+$v); print "$size\n";'`
	tail -n ${vt} all.tfeatures.all_x_vectors.${f} >> all.tfeatures.all_x_vectors.testing
	head -${v} all.tfeatures.all_y_vectors.${f} >> all.tfeatures.all_y_vectors.testing
done
