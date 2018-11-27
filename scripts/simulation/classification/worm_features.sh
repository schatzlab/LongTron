#!/bin/bash
#used to add sequence overlap features to a set of intervals (transcript, reads, etc...)
#this is primarily used for classification in terms of long read alignment error profiling
#features listed below should have already been preped separately (via prep_features.sh)

#Features (not in order):
#1) Splice motif count
#2) # of exons
#3) size of exons
#4) size of introns
#5) # of overlapping reads/transcripts
#6) # of common 150 SNPs overlapping

#samtools and bedtools should be in PATH
ST=samtools
BT=bedtools
PERBASE=/data3/LongReadRNA/scripts/simulation/classification/perbase
#the following needs to have had samtools faidx run on it ahead of time
GENOME_INDEX=c_elegans.PRJNA13758.WS267.genomic.fa
GENOME_SIZES=WS267.worm.chrom.sizes

SOURCE_PATH=./sources

GC=$SOURCE_PATH/gc.bg
RM=$SOURCE_PATH/rm.bed
SR=$SOURCE_PATH/sr.bed
SM=$SOURCE_PATH/splices.bed.bgz
EXONS_PERBASE=$SOURCE_PATH/exons.perbase.counts.bed.bgz
TRANSCRIPTS_PERBASE=$SOURCE_PATH/transcripts.perbase.counts.bed.bgz
LOCUS_STATS=$SOURCE_PATH/locus.stats.bed

#TODO?
#LOCAL_MAPPABILITY=$SOURCE_PATH/gv28.local_mappability.coords.bed
#GC=$SOURCE_PATH/gc5Base.bg.clean
#UMAP=$SOURCE_PATH/k24.Umap.MultiTrackMappability.sorted.bg
#SEGDUPS=$SOURCE_PATH/segmental_dups_hg38.sorted

BAM=$1

#convert long read alignment BAM into BED of features (# exons, exon length, intron length, min_exon, min_intron, mapping_quality)
#TODO: re-enable this for production
#samtools view -F 2308 $BAM | cut -f 1,2,3,4,5,6,10 | perl -ne 'BEGIN { open(IN,"<class2transcript.all"); while($line=<IN>) { chomp($line); ($class,$transcript)=split(/\t/,$line); $class{uc($transcript)}=$class; } close(IN); $b=0; open(IN,"<'${GENOME_SIZES}'"); while($line=<IN>) { ($c,$sz)=split(/\t/,$line); $chrms{$c}=1; } close(IN); } chomp; ($name,$flag,$c,$s,$mapping_quality,$f,$seq)=split(/\t/,$_); next if(!$chrms{$c}); $n=$name; $n=~s/^([^_]+).*$/$1/; $c1=$class{uc($n)}; $rl=length($seq); $start=$s-1; my $o = (int($flag) & 0x10)?"-":"+"; $nexons=1; $exon_length=0; $intron_length=0; $r=$start; $pr=$r; $min_exon_sz=2**31; $min_intron_sz=2**31; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N") { $nexons++; $el=$r-$pr; $il=$i; $min_exon_sz=$el if($el < $min_exon_sz); $min_intron_sz=$il if($il < $min_intron_sz); $exon_length+=$el; $intron_length+=$il; $pr=$r+$i; } $r+=$i; } $el=$r-$pr; $exon_length+=$el; $min_exon_sz=$el if($el < $min_exon_sz); $min_intron_sz=0 if($nexons == 1); print "".join("\t",($c,$start,$r,"r$b",0,$o,$n,$c1,$rl,$nexons,$exon_length,$intron_length,$min_exon_sz,$min_intron_sz,$mapping_quality))."\n"; $b++;' | sort -k1,1 -k2,2n -k3,3n > ${BAM}.bed.rl.nX3.minX2.mq

IN=${BAM}.bed.rl.nX3.minX2.mq

OFFSET=`head -1 ${IN} | tr \\\\t \\\\n | wc -l`

###RepeatMasker
#overlap with RepeatMasker, now we just get total overlap
$BT intersect -sorted -s -wao -a ${IN} -b $RM | perl -ne 'BEGIN { $bc=0; } chomp; @f=split(/\t/,$_); @f1=splice(@f,0,'$OFFSET'); $t=$f1[3]; $all=join("\t",@f1); if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $p=$all; } $pt=$t;  $rn=$f[3]; $rt=$f[6]; $bc+=$f[7]; END { if($pt) { print "$p\t$bc\n"; } }' > ${IN}.rm

###SimpleRepeats
$BT intersect -sorted -wao -a ${IN}.rm -b $SR | perl -ne 'BEGIN { $bc=0; } chomp; @f=split(/\t/,$_); @f1=splice(@f,0,('$OFFSET'+1)); $t=$f1[3]; $all=join("\t",@f1); if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $p=$all; } $pt=$t; $bc+=$f[3]; END { if($pt) { print "$p\t$bc\n"; } }' > ${IN}.rm.sr

IN=${IN}.rm.sr

###Count of overlapping transcripts/reads *on the same strand*, always has at least 1 (itself)
$BT intersect -sorted -s -c -a ${IN} -b ${IN} > ${IN}.ot

###exon density
cat ${IN}.ot | $PERBASE -c $GENOME_SIZES -f <(zcat $EXONS_PERBASE) > ${IN}.ot.ed

###transcript density
cat ${IN}.ot.ed | $PERBASE -c $GENOME_SIZES -f <(zcat $TRANSCRIPTS_PERBASE) > ${IN}.ot.ed.td

###Logs of nexons, exon bp, intron bp
cat ${IN}.ot.ed.td | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($ne,$ebp,$ibp)=($f[8],$f[9],$f[10]); $ibpl=($ibp>0?log($ibp):0); printf("%s\t%.3f\t%.3f\t%.3f\n",$f,log($ne),log($ebp),$ibpl);' > ${IN}.ot.ed.td.logsX3

###SpliceMotif frequency
#this is fixed now
cat ${IN}.ot.ed.td.logsX3 | $PERBASE -c $GENOME_SIZES -f <(zcat $SM) -m -s 5 >  ${IN}.ot.ed.td.logsX3.sm

###GC content
cat ${IN}.ot.ed.td.logsX3.sm | $PERBASE -c $GENOME_SIZES -f $GC -t c > ${IN}.ot.ed.td.logsX3.sm.gc

###Closest Locus transcript/exon stats (lengths, etc...), used at least for non-FL training/predictions
$BT closest -s -t first -D ref -a ${IN}.ot.ed.td.logsX3.sm.gc -b ${LOCUS_STATS} |  perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,('$OFFSET'+10)); $all=join("\t",@f1); @f2=splice(@f,6); $all2=join("\t",@f2); print "$all\t".$f[4]."\t$all2\n";' > ${IN}.ot.ed.td.logsX3.sm.gc.lsX10
