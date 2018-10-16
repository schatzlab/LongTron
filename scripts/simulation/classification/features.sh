#!/bin/bash
#used to add sequence overlap features to a set of intervals (transcript, reads, etc...)
#this is primarily used for classification in terms of long read alignment error profiling
#features listed below should have already been preped separately (via prep_features.sh)

#Features:
#1) Repeat Masker overlap
#2) Tandem Repeat Finder (TRF) overlap (subset of #1)
#3) Splice motif count
#4) # of exons
#5) size of exons
#6) size of introns
#7) # of overlapping reads/transcripts
#8) # of common 150 SNPs overlapping

#samtools and bedtools should be in PATH
ST=samtools
BT=bedtools
PERBASE=/data3/LongReadRNA/scripts/simulation/classification/perbase
#the following needs to have had samtools faidx run on it ahead of time
GENOME_INDEX=/data3/indexes/GRCh38_full_analysis_set_plus_decoy_hla.fa
GENOME_SIZES=/data/kent_tools/hg38.chrom.sizes.cannonical

SOURCE_PATH=/data7/schatz/simulation/sources

#RepeatMasker preped DB
RM=$SOURCE_PATH/rm.bed
#Simple Repeats preped DB
SR=$SOURCE_PATH/sr.bed
#Common SNPS 150 preped DB
SNPS=$SOURCE_PATH/snps.bed
#Splice Motifs
SM=$SOURCE_PATH/hg38_splice_motifs.all.bed.bgz
GC=$SOURCE_PATH/gc5Base.bg.clean
UMAP=$SOURCE_PATH/k24.Umap.MultiTrackMappability.sorted.bg
EXONS_PERBASE=$SOURCE_PATH/gencode.v28.basic.annotation.exons.perbase.counts.bgz
TRANSCRIPTS_PERBASE=$SOURCE_PATH/gencode.v28.basic.annotation.transcripts.perbase.counts.bgz

BAM=$1

#convert long read alignment BAM into BED of features (# exons, exon length, intron length, min_exon, min_intron, mapping_quality)
#TODO: re-enable this for production
samtools view -F 2308 $BAM | cut -f 1,2,3,4,5,6 | perl -ne 'BEGIN { $b=0; %chrms=("chr1"=>1,"chr2"=>1,"chr3"=>1,"chr4"=>1,"chr5"=>1,"chr6"=>1,"chr7"=>1,"chr8"=>1,"chr9"=>1,"chr10"=>1,"chr11"=>1,"chr12"=>1,"chr13"=>1,"chr14"=>1,"chr15"=>1,"chr16"=>1,"chr17"=>1,"chr18"=>1,"chr19"=>1,"chr20"=>1,"chr21"=>1,"chr22"=>1,"chrM"=>1,"chrX"=>1,"chrY"=>1); } chomp; ($name,$flag,$c,$s,$mapping_quality,$f)=split(/\t/,$_); next if(!$chrms{$c}); $start=$s--; my $o = (int($flag) & 0x10)?"-":"+"; $nexons=1; $exon_length=0; $intron_length=0; $r=$start; $pr=$r; $min_exon_sz=2**31; $min_intron_sz=2**31; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N") { $nexons++; $el=$r-$pr; $il=$i; $min_exon_sz=$el if($el < $min_exon_sz); $min_intron_sz=$il if($il < $min_intron_sz); $exon_length+=$el; $intron_length+=$il; $pr=$r+$i; } $r+=$i; } $el=$r-$pr; $exon_length+=$el; $min_exon_sz=$el if($el < $min_exon_sz); $min_intron_sz=0 if($nexons == 1); print "".join("\t",($c,$start,$r,"r$b.$name",0,$o,$nexons,$exon_length,$intron_length,$min_exon_sz,$min_intron_sz,$mapping_quality))."\n"; $b++;' | sort -k1,1 -k2,2n -k3,3n > ${BAM}.bed.n.min.mq

IN=${BAM}.bed.n.min.mq

OFFSET=12

###RepeatMasker
#overlap with RepeatMasker and then
#combine all the RMEs which match particular transcripts into one line counting total unique base overlap
#NOTE: we *do* track individual RMEs and their contribution to the overlap
$BT intersect -sorted -s -wao -a ${IN} -b $RM | perl -ne 'BEGIN { $bc=0; } chomp; @f=split(/\t/,$_); @f1=splice(@f,0,'$OFFSET'); $t=$f1[3]; $all=join("\t",@f1); if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $p=$all; } $pt=$t;  $rn=$f[3]; $rt=$f[6]; $bc+=$f[7]; END { if($pt) { print "$p\t$bc\n"; } }' > ${IN}.rm

###SimpleRepeats
$BT intersect -sorted -wao -a ${IN}.rm -b $SR | perl -ne 'BEGIN { $bc=0; } chomp; @f=split(/\t/,$_); @f1=splice(@f,0,('$OFFSET'+1)); $t=$f1[3]; $all=join("\t",@f1); if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $p=$all; } $pt=$t; $bc+=$f[3]; END { if($pt) { print "$p\t$bc\n"; } }' > ${IN}.rm.sr

###Common 150 SNPs
#NOTE: we don't track individual SNPs due to the potential size of the whole set
$BT intersect -sorted -s -wao -a ${IN}.rm.sr -b $SNPS | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,('$OFFSET'+2)); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $snps=""; $p=$all; } $pt=$t; if($f[6] != 0) { $bc++; $snps.=$f[3].";" } END { if($pt) { print "$p\t$bc\n"; } }' > ${IN}.rm.sr.snps

###Count of overlapping transcripts/reads *on the same strand*, always has at least 1 (itself)
$BT intersect -sorted -s -c -a ${IN}.rm.sr.snps -b ${IN} > ${IN}.rm.sr.snps.ot

###GC content
cat ${IN}.rm.sr.snps.ot | $PERBASE -c $GENOME_SIZES -f $GC -t c > ${IN}.rm.sr.snps.ot.gc

###Mappability (k=24, umap)
cat ${IN}.rm.sr.snps.ot.gc | $PERBASE -c $GENOME_SIZES -f $UMAP -t d > ${IN}.rm.sr.snps.ot.gc.umap


###Segmental Dups
#$BT intersect -sorted -wao -a ${IN}.rm.sr.snps.ot.gc.umap -b $SEGDUPS | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,('$OFFSET'+5)); ($s,$e)=($f1[1],$f1[2]); $d=($e-$s); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { $bc=$bc/$pd; printf("$p\t$c\t%.3f\n",$bc); } $bc=0; $c=0; $p=$all; } $pd=$d; $pt=$t; if($f[6] != 0) { $c++; $bc+=$f[6]; } END { if($pt) { $bc=$bc/$pd; printf("$p\t$c\t%.3f\n",$bc); } }' > ${IN}.rm.sr.snps.ot.gc.umap.nsd.sdo


###exon density
cat ${IN}.rm.sr.snps.ot.gc.umap | $PERBASE -c $GENOME_SIZES -f <(zcat $EXONS_PERBASE) > ${IN}.rm.sr.snps.ot.gc.umap.ed

###transcript density
cat ${IN}.rm.sr.snps.ot.gc.umap.ed | $PERBASE -c $GENOME_SIZES -f <(zcat $TRANSCRIPTS_PERBASE) > ${IN}.rm.sr.snps.ot.gc.umap.ed.td

###Logs of nexons, exon bp, intron bp
cat ${IN}.rm.sr.snps.ot.gc.umap.ed.td | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($ne,$ebp,$ibp)=($f[6],$f[7],$f[8]); $ibpl=($ibp>0?log($ibp):0); printf("%s\t%.3f\t%.3f\t%.3f\n",$f,log($ne),log($ebp),$ibpl);' > ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logs


###SpliceMotif frequency
#TODO: this is still a bit wonky (counts are off by a small(?, ~40) factor)
cat ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logs | $PERBASE -c $GENOME_SIZES -f <(zcat $SM) -t d -s 5 >  ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logs.sm


###need to fix the perbase version of this (counts are currently off by a small factor)
#still too slow way
#bedtools intersect -sorted -wao -s -a ${IN}.rm.sr.snps.ot -b <(zcat $SM) | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,13); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $p=$all; } $pt=$t; ($s1,$e1,$s2,$e2)=($f1[1],$f1[2],$f[1],$f[2]); $bc++ if($s2 >= $s1 && $e2 <= $e1); END { if($pt) { print "$p\t$bc\n"; }}' > ${IN}.rm.sr.snps.ot.sm
#old really slow way
#get splice motif frequency (note this was changed to use 1-base positions for querying via faidx)
#cat ${IN}.rm.sr.snps.ot | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($c,$s,$e,$n,$j,$o)=@f; $s++; @s=`'$ST' faidx '${GENOME_INDEX}' $c:$s-$e | fgrep -v ">"`; chomp(@s); $seq=join("",@s); @s=split(//,$seq); $nmotifs=0; $len=scalar(@s); $lm="GT"; $rm="AG"; if($o eq "-") { $lm="CT"; $rm="AC"; } for($i=0;$i+1 < $len;$i++) { $m=$s[$i].$s[$i+1]; if($m =~ /$lm/i || $m =~ /$rm/i) { $nmotifs++; } } print "$f\t$nmotifs\n";' > ${IN}.rm.sr.snps.ot.sm

