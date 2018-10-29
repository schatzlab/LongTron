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
SEGDUPS=$SOURCE_PATH/segmental_dups_hg38.sorted
EXONS_PERBASE=$SOURCE_PATH/gencode.v28.basic.annotation.exons.perbase.counts.bgz
TRANSCRIPTS_PERBASE=$SOURCE_PATH/gencode.v28.basic.annotation.transcripts.perbase.counts.bgz
LOCAL_MAPPABILITY=$SOURCE_PATH/gv28.local_mappability.coords.bed
LOCUS_STATS=$SOURCE_PATH/gencode.v28.basic.annotation.exons.stats.bed

BAM=$1

#convert long read alignment BAM into BED of features (# exons, exon length, intron length, min_exon, min_intron, mapping_quality)
#TODO: re-enable this for production
samtools view -F 2308 $BAM | cut -f 1,2,3,4,5,6,10 | perl -ne 'BEGIN { $b=0; %chrms=("chr1"=>1,"chr2"=>1,"chr3"=>1,"chr4"=>1,"chr5"=>1,"chr6"=>1,"chr7"=>1,"chr8"=>1,"chr9"=>1,"chr10"=>1,"chr11"=>1,"chr12"=>1,"chr13"=>1,"chr14"=>1,"chr15"=>1,"chr16"=>1,"chr17"=>1,"chr18"=>1,"chr19"=>1,"chr20"=>1,"chr21"=>1,"chr22"=>1,"chrM"=>1,"chrX"=>1,"chrY"=>1); } chomp; ($name,$flag,$c,$s,$mapping_quality,$f,$seq)=split(/\t/,$_); next if(!$chrms{$c}); $rl=length($seq); $start=$s-1; my $o = (int($flag) & 0x10)?"-":"+"; $nexons=1; $exon_length=0; $intron_length=0; $r=$start; $pr=$r; $min_exon_sz=2**31; $min_intron_sz=2**31; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N") { $nexons++; $el=$r-$pr; $il=$i; $min_exon_sz=$el if($el < $min_exon_sz); $min_intron_sz=$il if($il < $min_intron_sz); $exon_length+=$el; $intron_length+=$il; $pr=$r+$i; } $r+=$i; } $el=$r-$pr; $exon_length+=$el; $min_exon_sz=$el if($el < $min_exon_sz); $min_intron_sz=0 if($nexons == 1); print "".join("\t",($c,$start,$r,"r$b",0,$o,$name,$rl,$nexons,$exon_length,$intron_length,$min_exon_sz,$min_intron_sz,$mapping_quality))."\n"; $b++;' | sort -k1,1 -k2,2n -k3,3n > ${BAM}.bed.rl.nX3.minX2.mq

IN=${BAM}.bed.rl.nX3.minX2.mq

OFFSET=14

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


###exon density
cat ${IN}.rm.sr.snps.ot.gc.umap | $PERBASE -c $GENOME_SIZES -f <(zcat $EXONS_PERBASE) > ${IN}.rm.sr.snps.ot.gc.umap.ed

###transcript density
cat ${IN}.rm.sr.snps.ot.gc.umap.ed | $PERBASE -c $GENOME_SIZES -f <(zcat $TRANSCRIPTS_PERBASE) > ${IN}.rm.sr.snps.ot.gc.umap.ed.td

###Logs of nexons, exon bp, intron bp
cat ${IN}.rm.sr.snps.ot.gc.umap.ed.td | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($ne,$ebp,$ibp)=($f[7],$f[8],$f[9]); $ibpl=($ibp>0?log($ibp):0); printf("%s\t%.3f\t%.3f\t%.3f\n",$f,log($ne),log($ebp),$ibpl);' > ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3

###SpliceMotif frequency
#this is fixed now
cat ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3 | $PERBASE -c $GENOME_SIZES -f <(zcat $SM) -m -s 5 >  ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm

###Segmental Dups
$BT intersect -sorted -wao -a ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm -b $SEGDUPS | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,('$OFFSET'+12)); ($s,$e)=($f1[1],$f1[2]); $d=($e-$s); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { $bc=$bc/$pd; printf("$p\t$c\t%.3f\n",$bc); } $bc=0; $c=0; $p=$all; } $pd=$d; $pt=$t; if($f[6] != 0) { $c++; $bc+=$f[6]; } END { if($pt) { $bc=$bc/$pd; printf("$p\t$c\t%.3f\n",$bc); } }' > ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm.sdX2

###Local Mappability
#TODO: check offset on this one (OFFSET+14)
$BT intersect -sorted -s -wao -a ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm.sdX2 -b ${LOCAL_MAPPABILITY} | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,('$OFFSET'+14)); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { for($i=0;$i<scalar(@new);$i++) { $new[$i]/=$c; } $new_fields = join("\t",@new); print "$p\t$new_fields\n"; } $p=$all; @new=(); $c=0; } $c++; $pt=$t; for($i=6;$i<scalar(@f)-1;$i++) { $new[$i-6]+=$f[$i]; }  END { if($pt) { for($i=0;$i<scalar(@new);$i++) { $new[$i]/=$c; } $new_fields = join("\t",@new); print "$p\t$new_fields\n"; } }' > ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm.sdX2.lmX4

###Closest Locus transcript/exon stats (lengths, etc...), used at least for non-FL training/predictions
$BT closest -s -t first -D ref -a ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm.sdX2.lmX4 -b ${LOCUS_STATS} |  perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,('$OFFSET'+19)); $all=join("\t",@f1); @f2=splice(@f,6); $all2=join("\t",@f2); print "$all\t".$f[4]."\t$all2\n";' > ${IN}.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm.sdX2.lmX4.lsX9
