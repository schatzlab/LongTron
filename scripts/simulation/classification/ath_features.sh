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
GENOME_INDEX=TAIR10_Chr.all.fasta
GENOME_SIZES=TAIR10_Chr.all.fasta.chrm_sizes

SOURCE_PATH=./sources

SNPS=$SOURCE_PATH/snps.bed
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
samtools view -F 2308 $BAM | cut -f 1,2,3,4,5,6,10 | perl -ne 'BEGIN { open(IN,"<class2transcript.all"); while($line=<IN>) { chomp($line); ($class,$transcript)=split(/\t/,$line); $class{uc($transcript)}=$class; } close(IN); $b=0; %chrms=("Chr1"=>1,"Chr2"=>1,"Chr3"=>1,"Chr4"=>1,"Chr5"=>1,"ChrC"=>1,"ChrM"=>1); } chomp; ($name,$flag,$c,$s,$mapping_quality,$f,$seq)=split(/\t/,$_); next if(!$chrms{$c}); $n=$name; $n=~s/^([^\.]+\.\d+).*$/$1/; $c1=$class{uc($n)}; $rl=length($seq); $start=$s-1; my $o = (int($flag) & 0x10)?"-":"+"; $nexons=1; $exon_length=0; $intron_length=0; $r=$start; $pr=$r; $min_exon_sz=2**31; $min_intron_sz=2**31; while($f=~/(\d+)([NMD=X])/cg) { $i=$1; $t=$2; if($t eq "N") { $nexons++; $el=$r-$pr; $il=$i; $min_exon_sz=$el if($el < $min_exon_sz); $min_intron_sz=$il if($il < $min_intron_sz); $exon_length+=$el; $intron_length+=$il; $pr=$r+$i; } $r+=$i; } $el=$r-$pr; $exon_length+=$el; $min_exon_sz=$el if($el < $min_exon_sz); $min_intron_sz=0 if($nexons == 1); print "".join("\t",($c,$start,$r,"r$b",0,$o,$name,$c1,$rl,$nexons,$exon_length,$intron_length,$min_exon_sz,$min_intron_sz,$mapping_quality))."\n"; $b++;' | sort -k1,1 -k2,2n -k3,3n > ${BAM}.bed.rl.nX3.minX2.mq

IN=${BAM}.bed.rl.nX3.minX2.mq

OFFSET=`head -1 ${IN} | tr \\\\t \\\\n | wc -l`

###SNPs
#NOTE: we don't track individual SNPs due to the potential size of the whole set
$BT intersect -sorted -s -wao -a ${IN} -b $SNPS | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,('$OFFSET')); $all=join("\t",@f1); $t=$f1[3]; if($t ne $pt) { if($pt) { print "$p\t$bc\n"; } $bc=0; $snps=""; $p=$all; } $pt=$t; if($f[6] != 0) { $bc++; $snps.=$f[3].";" } END { if($pt) { print "$p\t$bc\n"; } }' > ${IN}.snps

###Count of overlapping transcripts/reads *on the same strand*, always has at least 1 (itself)
$BT intersect -sorted -s -c -a ${IN}.snps -b ${IN} > ${IN}.snps.ot

###exon density
cat ${IN}.snps.ot | $PERBASE -c $GENOME_SIZES -f <(zcat $EXONS_PERBASE) > ${IN}.snps.ot.ed

###transcript density
cat ${IN}.snps.ot.ed | $PERBASE -c $GENOME_SIZES -f <(zcat $TRANSCRIPTS_PERBASE) > ${IN}.snps.ot.ed.td

###Logs of nexons, exon bp, intron bp
cat ${IN}.snps.ot.ed.td | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($ne,$ebp,$ibp)=($f[8],$f[9],$f[10]); $ibpl=($ibp>0?log($ibp):0); printf("%s\t%.3f\t%.3f\t%.3f\n",$f,log($ne),log($ebp),$ibpl);' > ${IN}.snps.ot.ed.td.logsX3

###SpliceMotif frequency
#this is fixed now
cat ${IN}.snps.ot.ed.td.logsX3 | $PERBASE -c $GENOME_SIZES -f <(zcat $SM) -m -s 5 >  ${IN}.snps.ot.ed.td.logsX3.sm

###Closest Locus transcript/exon stats (lengths, etc...), used at least for non-FL training/predictions
$BT closest -s -t first -D ref -a ${IN}.snps.ot.ed.td.logsX3.sm -b ${LOCUS_STATS} |  perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); @f1=splice(@f,0,('$OFFSET'+8)); $all=join("\t",@f1); @f2=splice(@f,6); $all2=join("\t",@f2); print "$all\t".$f[4]."\t$all2\n";' > ${IN}.snps.ot.ed.td.logsX3.sm.lsX10
