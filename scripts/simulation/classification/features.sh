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
#the following needs to have had samtools faidx run on it ahead of time
GENOME_INDEX=/data3/indexes/GRCh38_full_analysis_set_plus_decoy_hla.fa
GENOME_SIZES=/data/kent_tools/hg38.chrom.sizes.cannonical

#RepeatMasker preped DB
RM=rm.bed
#Simple Repeats preped DB
SR=sr.bed
#Common SNPS 150 preped DB
SNPS=snps.bed
#Splice Motifs
SM=hg38_splice_motifs.all.bed.bgz
GC=gc5Base.bg.clean
UMAP=k24.Umap.MultiTrackMappability.sorted.bg
EXONS_PERBASE=gencode.v28.basic.annotation.exons.perbase.counts.bgz
TRANSCRIPTS_PERBASE=gencode.v28.basic.annotation.transcripts.perbase.counts.bgz

#assume we get a true BED file as input
IN=$1

###Exon, intron stats
#BED format: 
#chr1    65507   71582   r0      0       +       65507   71582   255,0,0 2       66,2546 0,3529
cat $IN | perl -ne 'chomp; $f=$_; ($c,$s,$e,$n,$j,$o,$ts,$te,$b,$nb,$szs,$starts)=split(/\t/,$f); $nexons=$nb; $nintrons=$nexons-1; $exon_bp=0; map { $exon_bp+=$_; } split(/,/,$szs); $total_size = ($e-$s); $intron_bp = $total_size - $exon_bp; print "".(join("\t",($c,$s,$e,$n,0,$o,$nexons,$exon_bp,$intron_bp)))."\n";' | sort -k1,1 -k2,2n -k3,3n > ${IN}.n

###RepeatMasker
#overlap with RepeatMasker and then
#combine all the RMEs which match particular transcripts into one line counting total unique base overlap
#NOTE: we *do* track individual RMEs and their contribution to the overlap
$BT intersect -sorted -s -wao -a ${IN}.n -b $RM | perl -ne 'BEGIN { $bc=0; } chomp; @f=split(/\t/,$_); @f1=splice(@f,0,9); $t=$f1[3]; $all=join("\t",@f1); if($t ne $pt) { if($pt) { $p=~s/\t;/\t/; print "$p\t$bc\n"; } $bc=0; $p=$all."\t"; } $pt=$t;  $rn=$f[3]; $rt=$f[6]; $p.=";$rn:$rt:".$f[7]; $bc+=$f[7]; END { if($pt) { $p=~s/\t;/\t/; print "$p\t$bc\n"; } }' > ${IN}.n.rm

###SimpleRepeats
$BT intersect -sorted -wao -a ${IN}.n.rm -b $SR | perl -ne 'BEGIN { $bc=0; } chomp; @f=split(/\t/,$_); @f1=splice(@f,0,11); $t=$f1[3]; $all=join("\t",@f1); if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $p=$all; } $pt=$t; $bc+=$f[3]; END { if($pt) { print "$p\t$bc\n"; } }' > ${IN}.n.rm.sr

###Common 150 SNPs
#NOTE: we don't track individual SNPs due to the potential size of the whole set
$BT intersect -sorted -s -wao -a ${IN}.n.rm.sr -b $SNPS | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,12); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $snps=""; $p=$all; } $pt=$t; if($f[6] != 0) { $bc++; $snps.=$f[3].";" } END { if($pt) { print "$p\t$bc\n"; } }' > ${IN}.n.rm.sr.snps

###Count of overlapping transcripts/reads *on the same strand*, always has at least 1 (itself)
$BT intersect -sorted -s -c -a ${IN}.n.rm.sr -b ${IN}.n > ${IN}.n.rm.sr.snps.ot

###GC content
cat ${IN}.n.rm.sr.snps.ot | $PERBASE -c $GENOME_SIZES -f $GC -t c > ${IN}.n.rm.sr.snps.ot.gc

###Mappability (k=24, umap)
cat ${IN}.n.rm.sr.snps.ot.gc | $PERBASE -c $GENOME_SIZES -f $UMAP -t d > ${IN}.n.rm.sr.snps.ot.gc.umap


###Segmental Dups
$BT intersect -sorted -wao -a ${IN}.n.rm.sr.snps.ot.gc.umap -b $SEGDUPS | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,18); ($s,$e)=($f1[1],$f1[2]); $d=($e-$s); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { $bc=$bc/$pd; printf("$p\t$c\t%.3f\n",$bc); } $bc=0; $c=0; $p=$all; } $pd=$d; $pt=$t; if($f[6] != 0) { $c++; $bc+=$f[6]; } END { if($pt) { $bc=$bc/$pd; printf("$p\t$c\t%.3f\n",$bc); } }' > ${IN}.n.rm.sr.snps.ot.gc.umap.nsd.sdo


###exon density
cat ${IN}.n.rm.sr.snps.ot.gc.umap.nsd.sdo | $PERBASE -c $GENOME_SIZES -f <(zcat $EXONS_PERBASE) > ${IN}.n.rm.sr.snps.ot.gc.umap.nsd.sdo.ed

###transcript density
cat ${IN}.n.rm.sr.snps.ot.gc.umap.nsd.sdo.ed | $PERBASE -c $GENOME_SIZES -f <(zcat $TRANSCRIPTS_PERBASE) > ${IN}.n.rm.sr.snps.ot.gc.umap.nsd.sdo.ed.td

###Logs of nexons, exon bp, intron bp
cat ${IN}.n.rm.sr.snps.ot.gc.umap.nsd.sdo.ed.td | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($ne,$ebp,$ibp)=($f[8],$f[9],$f[10]); $idbpl=($ibp>0?log($ibp):0); printf("%s\t%.3f\t%.3f\t%.3f\n",$f,log($ne),log($ebp),$ibpl);' > ${IN}.n.rm.sr.snps.ot.gc.umap.nsd.sdo.ed.td.logs


###SpliceMotif frequency
###need to fix the perbase version of this (counts are currently off by a small factor)
#still too slow way
bedtools intersect -sorted -wao -s -a ${IN}.n.rm.sr.snps.ot -b <(zcat $SM) | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,13); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $p=$all; } $pt=$t; ($s1,$e1,$s2,$e2)=($f1[1],$f1[2],$f[1],$f[2]); $bc++ if($s2 >= $s1 && $e2 <= $e1); END { if($pt) { print "$p\t$bc\n"; }}' > ${IN}.n.rm.sr.snps.ot.sm
#old really slow way
#get splice motif frequency (note this was changed to use 1-base positions for querying via faidx)
cat ${IN}.n.rm.sr.snps.ot | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($c,$s,$e,$n,$j,$o)=@f; $s++; @s=`'$ST' faidx '${GENOME_INDEX}' $c:$s-$e | fgrep -v ">"`; chomp(@s); $seq=join("",@s); @s=split(//,$seq); $nmotifs=0; $len=scalar(@s); $lm="GT"; $rm="AG"; if($o eq "-") { $lm="CT"; $rm="AC"; } for($i=0;$i+1 < $len;$i++) { $m=$s[$i].$s[$i+1]; if($m =~ /$lm/i || $m =~ /$rm/i) { $nmotifs++; } } print "$f\t$nmotifs\n";' > ${IN}.n.rm.sr.snps.ot.sm

