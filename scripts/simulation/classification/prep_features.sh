#!/bin/bash
#used to add sequence overlap features to a set of intervals (transcript, reads, etc...)
#this is primarily used for classification in terms of long read alignment error profiling
#assumes bedtools and samtools are in the PATH

#Features:
#1) Repeat Masker overlap
#2) Tandem Repeat Finder (TRF) overlap (subset of #1)
#3) Splice motif count
#4) # of exons
#5) size of exons
#6) size of introns
#7) # of overlapping reads/transcripts
#8) # of common 150 SNPs overlapping

ANNOT_VER=HG38
CHROM_SIZES=/data/kent_tools/hg38.chrom.sizes
#number of cannonical chromosomes (e.g. in human 1-22, X,Y,M)
NUM_CHRMS=25
ANNOTATION=gencode.v28.basic.annotation.gtf.gz
SNP=snp150Common.txt.gz
RM=hg38_repeatmasker_rmsk.gz
GC=gc5Base.bw
MAPPABILITY=k24.Umap.MultiTrackMappability.bw

BW2BG=/data/kent_tools/bigWigToBedGraph

#pull out basic per-transcript/exon stats from annotation
#this does NOT produce a BED file, coordinates are still 1-base
zcat ${ANNOTATION} | egrep -e '	exon	' | sort -k1,1 -k4,4n -k5,5n > ${ANNOTATION}.exons
cat ${ANNOTATION}.exons | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($c,$s,$e,$o)=($f[0],$f[3],$f[4],$f[6]); $f=~/gene_id "([^"]+)";/; $g=$1; $f=~/gene_type "([^"]+)";/; $gt=$1; $f=~/gene_name "([^"]+)";/; $gn=$1; $f=~/transcript_id "([^"]+)";/; $t=$1; $ts{$t}->{gene}=[$g,$gn,$gt,$c,$o]; $d=($e-$s)+1; $ts{$t}->{esum}+=$d; push(@{$ts{$t}->{exons}},[$s,$e,$d]); END { for $t (keys %ts) { $v=$ts{$t}; ($g,$gn,$gt,$c,$o)=@{$v->{gene}}; $esum=$v->{esum}; $ne=scalar(@{$v->{exons}}); $ni=0; $isum=0; $exons=""; $introns=""; $tstart=2**31; $tend=-1; $pe=undef; $min_exon_sz=2**31; $min_intron_sz=2**31; for $exon (@{$v->{exons}}) {  ($s,$e)=@$exon; $sz=($e-$s)+1; $min_exon_sz=$sz if($sz < $min_exon_sz); $tstart=$s if($s < $tstart); $tend=$e if($e > $tend); $exons.="$s-$e;"; if($pe) { $is=$pe+1; $ie=$s-1; $ni++; $sz=($ie-$is)+1; $min_intron_sz=$sz if($sz < $min_intron_sz); $isum+=($ie-$is)+1; $introns.="$is-$ie;"; } $pe=$e; } $exons=~s/;$//; $introns=~s/;$//; $min_intron_sz=0 if($ne==1); print "".(join("\t",($c,$tstart,$tend,$o,$t,$g,$gn,$gt,$ne,$esum,$ni,$isum,$exons,$introns,$min_exon_sz,$min_intron_sz)))."\n"; }}' > ${ANNOTATION}.exons_introns_per_transcript

#pull out annotated junctions
cat ${ANNOTATION}.exons | cut -f 1,4,5,7,9 | perl -ne 'chomp; $f=$_; ($c,$s,$e,$o,$info)=split(/\t/,$f); $info=~/transcript_id "([^"]+)";/; $t=$1; $info=~/exon_number (\d+);/; $eid=$1; $t.= "_0_$o"; if($ts{$t}) { $pe=$ts{$t}; $c1=$pe+1; $c2=$s-1;  $eid--; print "$c\t$c1\t$c2\t$t\t$eid\n"; } $ts{$t}=$e;' | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; $f=$_; ($c,$s,$e,$n)=split(/\t/,$f); if($n=~/_\+$/) { print "$f\n"; next; } if(!$ns{$n}) { $eid=0; } else { $eid=$ns{$n} } $eid++; print "$c\t$s\t$e\t$n\t$eid\n"; $ns{$n}=$eid;'  > ${ANNOTATION}.junctions

#get number of exon per transcript
cat ${ANNOTATION}.exons | perl -ne 'chomp; $f=$_; $f=~/transcript_id "([^"]+)";/; $t=$1; print "$t\n";' | sort | uniq -c | sort -k1,1nr | perl -ne 'chomp; ($j,$c,$n)=split(/\s+/,$_); print "$n\t$c\n";' > ${ANNOTATION}.transcripts_exon_count


#convert SNPs to BED
zcat $SNPS | cut -f 2-5,7 | sort -k1,1 -k2,2n -k3,3n > ${SNPS}.sorted

#Repeat Makser & SNPs prep
#sort, collapse overlaps by strand, and then sort the combined file again
zcat ${RM} | sort -k1,1 -k2,2n -k3,3n > ${RM}.sorted
echo "" > ${RM}.strand
echo "" > ${SNPS}.strand
for sign in '\+' '-';
do
	egrep -e "	${sign}	" ${RM}.sorted | perl -ne 'chomp; $f=$_; ($c,$s,$e,$o,$n,$t)=split(/\t/,$f); if($p && $pc eq $c && $s <= $pe) { $p.=";$n"; $pt.=";$t"; $pe=$e; next; } if($p) { print "$pc\t$ps\t$pe\t$po\t$p\t$pt\n";} $p=$n; $pt=$t; $pc=$c; $ps=$s; $pe=$e; $po=$o; END { if($p) { print "$pc\t$ps\t$pe\t$po\t$p\t$pt\n";  } }' >> ${RM}.strand
	egrep -e "	${sign}$" ${SNPS}.sorted | perl -ne 'chomp; $f=$_; ($c,$s,$e,$n,$o)=split(/\t/,$f); if($p && $pc eq $c && $s+1 <= $pe) { $p.=";$n"; $pe=$e; next; } if($p) { print "$pc\t$ps\t$pe\t$p\t$po\n";} $p=$n; $pc=$c; $ps=$s; $pe=$e; $po=$o; END { if($p) { print "$pc\t$ps\t$pe\t$p\t$po\n";  } }' >> ${SNPS}.strand
done
cat ${RM}.strand | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; ($c,$s,$e,$n,$j,$o,$t)=split(/\t/,$_); $s--; print "".join("\t",($c,$s,$e,$n,$j,$o,$t))."\n";' > rm.bed
#TODO: get rid of "bad" chromosomes in SNPs
cat ${SNPS}.strand | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; $f=$_; ($c,$s,$e,$n,$o)=split(/\t/,$f); $s--; print "".join("\t",($c,$s,$e,$n,0,$o))."\n";' > snps.bed

#simple (TRF) repeats
zcat simple_repeats_hg38.bed.gz | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; $f=$_; ($c,$s,$e)=split(/\t/,$f); if($pc && $pc eq $c && $s+1 <= $pe) { $pc=$c; $pe=$e; next; } if($pc) { print "$pc\t$ps\t$pe\n"; } $pc=$c; $ps=$s; $pe=$e; END { if($pc) { print "$pc\t$ps\t$pe\n";  } }' > sr.bed

#get mapping quality per read alignment
samtools view -F 2308 trans_sim10.fl.bam | cut -f 1,5 | perl -ne 'BEGIN { $b=0; } chomp; ($t,$mq)=split(/\t/,$_); @f=split(/_/,$t); $t=shift(@f); $t=~s/\.[+-]$//; print "r$b\t$t\t$mq\n"; $b++;' > trans_sim10.fl.bam.bed2.mapqualitys

#prep GC/UMAP track data by removing non-cannonical chromosomes and sort
#assume we're getting BedGraph files (bigWigToBedGraph has already been run)
$BW2BG $GC ${GC}.bg
cat ${GC}.bg | perl -ne 'chomp; $f=$_; @f=split(/\t/,$_); next if($f[0]=~/(_random)|(_alt)|(chrUn)/i); print "$f\n";' | sort -k1,1 -k2,2n -k3,3n > ${GC}.clean.sorted.bg
$BW2BG ${MAPPABILITY} ${MAPPABILITY}.bg
cat ${MAPPABILITY}.bg | perl -ne 'chomp; $f=$_; @f=split(/\t/,$_); next if($f[0]=~/(_random)|(_alt)|(chrUn)/i); print "$f\n";' | sort -k1,1 -k2,2n -k3,3n > ${MAPPABILITY}.clean.sorted.bg


#create perbase for exon/transcript densities
#these are big files, so have the ANNOTATION file sitting on a large filesystem with plenty of free space
zcat $ANNOTATION | egrep '	exon	' | cut -f 1,4,5,7 | sort -k1,1 -k2,2n -k3,3n | perl -ne 'BEGIN { open(OUT,">'${ANNOTATION}.exons.sorted.bed'"); } chomp; $f=$_; ($c,$s,$e,$o)=split(/\t/,$f); print OUT "$c\t".($s-1)."\t$e\n"; if(!$pc || ($pc ne $c || $s > $pe)) { if($pc) { $ps--; for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; }} $ps=$s; $pc=$c; $pe=$e; $po=$o; next; } $pe=$e; END { if($pc) { $ps--; for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; } } close(OUT);}' | bgzip > ${$ANNOTATION}.exons.perbase.bgz

bedtools intersect -sorted -c -a <(zcat ${$ANNOTATION}.exons.perbase.bgz) -b <(sort -k1,1 -k2,2n -k3,3n ${ANNOTATION}.exons.sorted.bed) | bgzip > ${ANNOTATION}.exons.perbase.counts.bgz

zcat $ANNOTATION | egrep '	transcript	' | cut -f 1,4,5,7 | sort -k1,1 -k2,2n -k3,3n | perl -ne 'BEGIN { open(OUT,">'${ANNOTATION}.transcripts.sorted.bed'"); } chomp; $f=$_; ($c,$s,$e,$o)=split(/\t/,$f); print OUT "$c\t".($s-1)."\t$e\n"; if(!$pc || ($pc ne $c || $s > $pe)) { if($pc) { $ps--; for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; }} $ps=$s; $pc=$c; $pe=$e; $po=$o; next; } $pe=$e; END { if($pc) { $ps--; for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; } } close(OUT);}' | bgzip ${$ANNOTATION}.transcripts.perbase.bgz

bedtools intersect -sorted -c -a <(zcat ${$ANNOTATION}.transcripts.perbase.bgz) -b <(sort -k1,1 -k2,2n -k3,3n ${ANNOTATION}.transcripts.sorted.bed) | bgzip > ${ANNOTATION}.transcripts.perbase.counts.bgz

#extract out all cannonical splice motifs (forward + reverese) from reference across all main chromosomes
#takes ~33m22.814s for HG38
(head -${NUM_CHRMS} $CHROM_SIZES | cut -f 1 | sort | perl -ne 'BEGIN { %F=("GT"=>1,"gt"=>1,"Gt"=>1,"gT"=>1,"AG"=>1,"ag"=>1,"Ag"=>1,"aG"=>1); %R=("CT"=>1,"AC"=>1,"cT"=>1,"Ct"=>1,"ct"=>1,"Ac"=>1,"aC"=>1,"ac"=>1); } chomp; $chrm=$_; @s=`samtools faidx /data3/indexes/GRCh38_full_analysis_set_plus_decoy_hla.fa $chrm | fgrep -v ">"`; chomp(@s); $seq=join("",@s); $len=length($seq); for($i=0;$i<($len-1);$i++) { $m=substr($seq,$i,2); $b="$chrm\t$i\t".($i+1)."\t$m\n"; if($F{$m}) { print $b; } if($R{$m}) { print STDERR $b; } }' | bgzip > ${ANNOT_VER}.forward_splice_motifs.all.tsv.bgz) 2>&1 | bgzip > ${ANNOT_VER}.reverse_splice_motifs.all.tsv.bgz

#merge forward and reverse strands of all cannonical splice motifs from reference into one file (sorted)
cat <(zcat ${ANNOT_VER}.forward_splice_motifs.all.tsv.bgz | perl -ne 'chomp; print "$_\t0\t+\n";') <(zcat ${ANNOT_VER}.reverse_splice_motifs.all.tsv.bgz | perl -ne 'chomp; print "$_\t0\t-\n";') | sort -k1,1 -k2,2n -k3,3n | bgzip > ${ANNOT_VER}.splice_motifs.all.tsv.bgz
