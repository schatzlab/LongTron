#!/bin/bash
#generate a set of 1) training 2) validation and 3) test partitions for each category of transcript/reads

ntraining=$1
nvalidation=$2

#prep transcript features
#name,nexons,exon_bp,intron_bp,rm_overlap,simpler_overlap,common_snp150_count,splice_motif_count,nTranscripts,nTSS,nPolyA,gene_smotifs
#nTranscripts,nTSS,nPolyA have little feature importance
#cut -f 1,11,12,14,18,19,21- all.tfeatures | tail -n+2 | perl -ne 'BEGIN { %h=("non-recurrent"=>0,"novel"=>1,"problem-free"=>2,"recurrent"=>3); } chomp; @f=split(/\t/,$_); $n=shift(@f); @a=(0,0,0,0); $a[$h{$n}]=1; print "$n\t".join("\t",@f)."\n"; print STDERR "$n\t".join("\t",@a)."\n";' > features.x 2> features.y

#prep *read alignments* features
#cut -f4,7-9,11- trans_sim10.fl.bam.bed.n.rm.sr.ot.sm.snps | perl -ne 'BEGIN { %c=("non-recurrent"=>0,"novel"=>1,"problem-free"=>2,"recurrent"=>3); open(IN,"<transcript2category2rid.mapping"); while($line=<IN>) { chomp($line); ($t,$c,$r)=split(/\s+/,$line); $h{$r}=$c; } close(IN); } chomp; @f=split(/\t/,$_); $n=shift(@f); $c=$h{$n}; @a=(0,0,0,0); $a[$c{$c}]=1; print "$c\t".join("\t",@f)."\n"; print STDERR "$c\t".join("\t",@a)."\n";' > trans_sim10.fl.bam.bed.n.rm.sr.ot.sm.snps.features 2>trans_sim10.fl.bam.bed.n.rm.sr.ot.sm.snps.labels

ln -s trans_sim10.fl.bam.bed.n.rm.sr.ot.sm.snps.features features.x
ln -s trans_sim10.fl.bam.bed.n.rm.sr.ot.sm.snps.labels features.y

echo -n "" > training.x
echo -n "" > training.y
echo -n "" > validation.x
echo -n "" > validation.y
echo -n "" > testing.x
echo -n "" > testing.y
echo -n "" > num_training_samples_per_cat
for f in "non-recurrent" "novel" "problem-free" "recurrent";
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
for f in "non-recurrent" "novel" "problem-free" "recurrent";
do
	head -${min} features.x.${f} >> training.x.bal
	head -${min} features.y.${f} >> training.y.bal
done


#create 1) binary (0,1) labels and separately 2) only use recurrent and problem-free
n=`head -1 training.x | tr \\\t \\\n | wc -l | cut -d" " -f 1`
for t in "training" "validation" "testing";
do
	cat ${t}.y | perl -ne 'chomp; ($non_recurrent,$novel,$problem_free,$recurrent) = split(/\t/,$_); if($problem_free == 1) { print "0\n"; next; } print "1\n";' > binary.y.${t}
	#binary: recurrent vs. problem-free
	paste ${t}.x ${t}.y | egrep -e '	0	0	.	.$' | perl -ne 'chomp; @f=split(/\t/,$_); @f1=splice(@f,0,'${n}'); print "".join("\t",@f1)."\n"; ($non_recurrent,$novel,$problem_free,$recurrent) = @f; print STDERR "$recurrent\n";' > ${t}.2.x 2> ${t}.2.y
	#binary: recurrent + non-recurrent vs. problem-free
	paste ${t}.x ${t}.y | egrep -e '	.	0	.	.$' | perl -ne 'chomp; @f=split(/\t/,$_); @f1=splice(@f,0,'${n}'); print "".join("\t",@f1)."\n"; ($non_recurrent,$novel,$problem_free,$recurrent) = @f; print STDERR "$problem_free\n";' > ${t}.3a.x 2> ${t}.3a.y
	#binary: recurrent + novels vs. problem-free
	paste ${t}.x ${t}.y | egrep -e '	0	.	.	.$' | perl -ne 'chomp; @f=split(/\t/,$_); @f1=splice(@f,0,'${n}'); print "".join("\t",@f1)."\n"; ($non_recurrent,$novel,$problem_free,$recurrent) = @f; print STDERR "$problem_free\n";' > ${t}.3b.x 2> ${t}.3b.y
done
