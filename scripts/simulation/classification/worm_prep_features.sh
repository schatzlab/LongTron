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

ANNOT_VER=WS267
CHROM_SIZES=${ANNOT_VER}.worm.chrom.sizes
#number of cannonical chromosomes (e.g. in human 1-22, X,Y,M)
NUM_CHRMS=7
ANNOTATION=c_elegans.PRJNA13758.WS267.canonical_geneset.gtf.gz
GENOME=c_elegans.PRJNA13758.WS267.genomic.fa
#SNP=snp150Common.txt.gz
RM=rmsk.txt.gz
SR=simpleRepeat.txt.gz
GC=gc5Base.bw
#MAPPABILITY=k24.Umap.MultiTrackMappability.bw

BW2BG=/data/kent_tools/bigWigToBedGraph

#get chrm sizes
cat $GENOME | perl -ne 'chomp; $s=$_; if($s=~/^>/) { if($seq && length($seq) > 0) { print "\t".length($seq)."\n"; } $seq=""; $s=~s/^>//; print "$s"; next; } $seq.=$s; END { if($seq && length($seq) > 0) { print "\t".length($seq)."\n"; } }' > $CHROM_SIZES

NUM_CHRMS=`wc -l $CHROM_SIZES | cut -d' ' -f 1`

#pull out basic per-transcript/exon stats from annotation
#this does NOT produce a BED file, coordinates are still 1-base
zcat ${ANNOTATION} | egrep -e '	exon	' | sort -k1,1 -k4,4n -k5,5n > ${ANNOTATION}.exons

#pull out annotated junctions
cat ${ANNOTATION}.exons | cut -f 1,4,5,7,9 | perl -ne 'BEGIN { open(OUT,"| sort -k1,1 -k2,2n -k3,3n >dmel-all-r6.23.gtf.exons.bed"); } chomp; $f=$_; ($c,$s,$e,$o,$info)=split(/\t/,$f); $info=~/transcript_id "([^"]+)";/; $t=$1; print OUT "$c\t".($s-1)."\t$e\t$t\t0\t$o\n"; if($ts{$t}) { $pe=$ts{$t}; $c1=$pe+1; $c2=$s-1;  print "$c\t$c1\t$c2\t$t\n"; } $ts{$t}=$e;' | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; $f=$_; ($c,$s,$e,$n)=split(/\t/,$f); if(!$ns{$n}) { $eid=0; } else { $eid=$ns{$n} } $eid++; print "$c\t$s\t$e\t$n\t$eid\n"; $ns{$n}=$eid; END { close(IN); }'  >  ${ANNOTATION}.junctions

#get number of exon per transcript
cat ${ANNOTATION}.exons | perl -ne 'chomp; $f=$_; $f=~/transcript_id "([^"]+)";/; $t=$1; print "$t\n";' | sort | uniq -c | sort -k1,1nr | perl -ne 'chomp; ($j,$c,$n)=split(/\s+/,$_); print "".uc($n)."\t$c\n";' > ${ANNOTATION}.transcripts_exon_count

#Repeat Makser & SNPs prep
#sort, collapse overlaps by strand, and then sort the combined file again
zcat ${RM} | cut -f 6-8,10-12 | sort -k1,1 -k2,2n -k3,3n > ${RM}.sorted
for sign in '\+' '-';
do
	egrep -e "	${sign}	" ${RM}.sorted | perl -ne 'chomp; $f=$_; ($c,$s,$e,$o,$n,$t)=split(/\t/,$f); $c=~/^chr//; if($p && $pc eq $c && $s <= $pe) { $p.=";$n"; $pt.=";$t"; $pe=$e; next; } if($p) { print "$pc\t$ps\t$pe\t$po\t$p\t$pt\n";} $p=$n; $pt=$t; $pc=$c; $ps=$s; $pe=$e; $po=$o; END { if($p) { print "$pc\t$ps\t$pe\t$po\t$p\t$pt\n";  } }' > ${RM}.strand
done
cat ${RM}.strand | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; @f=split(/\t/,$_); ($c,$s,$e)=splice(@f,0,3); $s--; print "".join("\t",($c,$s,$e))."\t".join("\t",@f)."\n";' > rm.bed

#simple (TRF) repeats
zcat $SR | cut -f 2- | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; $f=$_; ($c,$s,$e)=split(/\t/,$f); $c=~/^chr//; if($pc && $pc eq $c && $s+1 <= $pe) { $pc=$c; $pe=$e; next; } if($pc) { print "$pc\t$ps\t$pe\n"; } $pc=$c; $ps=$s; $pe=$e; END { if($pc) { print "$pc\t$ps\t$pe\n";  } }' > sr.bed

#prep GC/UMAP track data by removing non-cannonical chromosomes and sort
#assume we're getting BedGraph files (bigWigToBedGraph has already been run)
$BW2BG $GC ${GC}.bg
#cat ${GC}.bg | perl -ne 'chomp; $f=$_; @f=split(/\t/,$_); next if($f[0]=~/(_random)|(_alt)|(chrUn)/i); print "$f\n";' | sort -k1,1 -k2,2n -k3,3n > ${GC}.clean.sorted.bg
mv ${GC}.bg ${GC}.clean.sorted.bg

#create perbase for exon/transcript densities
#this assumes that there are multiple overlapping exons in some regions (maybe even duplicates across transcripts)
#these are big files, so have the ANNOTATION file sitting on a large filesystem with plenty of free space
zcat $ANNOTATION | egrep '	exon	' | cut -f 1,4,5,7 | sort -k1,1 -k2,2n -k3,3n | perl -ne 'BEGIN { open(OUT,">'${ANNOTATION}'.exons.sorted.bed"); } chomp; $f=$_; ($c,$s,$e,$o)=split(/\t/,$f); print OUT "$c\t".($s-1)."\t$e\n"; if(!$pc || ($pc ne $c || $s > $pe)) { if($pc) { $ps--; for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; }} $ps=$s; $pc=$c; $pe=$e; $po=$o; next; } $pe=$e; END { if($pc) { $ps--; for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; } } close(OUT);}' | bgzip > ${ANNOTATION}.exons.perbase.bgz

bedtools intersect -sorted -c -a <(zcat ${ANNOTATION}.exons.perbase.bgz) -b <(sort -k1,1 -k2,2n -k3,3n ${ANNOTATION}.exons.sorted.bed) | bgzip > ${ANNOTATION}.exons.perbase.counts.bgz

zcat $ANNOTATION | egrep '	transcript	' | cut -f 1,4,5,7 | sort -k1,1 -k2,2n -k3,3n | perl -ne 'BEGIN { open(OUT,">'${ANNOTATION}'.transcripts.sorted.bed"); } chomp; $f=$_; ($c,$s,$e,$o)=split(/\t/,$f); print OUT "$c\t".($s-1)."\t$e\n"; if(!$pc || ($pc ne $c || $s > $pe)) { if($pc) { $ps--; for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; }} $ps=$s; $pc=$c; $pe=$e; $po=$o; next; } $pe=$e; END { if($pc) { $ps--; for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; } } close(OUT);}' | bgzip > ${ANNOTATION}.transcripts.perbase.bgz

bedtools intersect -sorted -c -a <(zcat ${ANNOTATION}.transcripts.perbase.bgz) -b <(sort -k1,1 -k2,2n -k3,3n ${ANNOTATION}.transcripts.sorted.bed) | bgzip > ${ANNOTATION}.transcripts.perbase.counts.bgz

#extract out all cannonical splice motifs (forward + reverese) from reference across all main chromosomes
#takes ~33m22.814s for HG38
(head -${NUM_CHRMS} $CHROM_SIZES | cut -f 1 | sort | perl -ne 'BEGIN { %F=("GT"=>1,"gt"=>1,"Gt"=>1,"gT"=>1,"AG"=>1,"ag"=>1,"Ag"=>1,"aG"=>1); %R=("CT"=>1,"AC"=>1,"cT"=>1,"Ct"=>1,"ct"=>1,"Ac"=>1,"aC"=>1,"ac"=>1); } chomp; $chrm=$_; @s=`samtools faidx '${GENOME}' $chrm | fgrep -v ">"`; chomp(@s); $seq=join("",@s); $len=length($seq); for($i=0;$i<($len-1);$i++) { $m=substr($seq,$i,2); $b="$chrm\t$i\t".($i+1)."\t$m\n"; if($F{$m}) { print $b; } if($R{$m}) { print STDERR $b; } }' | bgzip > ${ANNOT_VER}.forward_splice_motifs.all.tsv.bgz) 2>&1 | bgzip > ${ANNOT_VER}.reverse_splice_motifs.all.tsv.bgz

#merge forward and reverse strands of all cannonical splice motifs from reference into one file (sorted)
cat <(zcat ${ANNOT_VER}.forward_splice_motifs.all.tsv.bgz | perl -ne 'chomp; print "$_\t0\t+\n";') <(zcat ${ANNOT_VER}.reverse_splice_motifs.all.tsv.bgz | perl -ne 'chomp; print "$_\t0\t-\n";') | sort -k1,1 -k2,2n -k3,3n | bgzip > ${ANNOT_VER}.splice_motifs.all.tsv.bgz


#get per-gene level stats about transcript and exon lengths/counts/sums/averages
zcat $ANNOTATION | egrep -e '	exon	' | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($ch,$start,$end,$strand)=($f[0],$f[3],$f[4],$f[6]); $f=~/gene_id "([^"]+)"/; $g=$1; $h2{$g}->[0]=$start if(!$h2{$g}->[0] || $start < $h2{$g}->[0]); $h2{$g}->[1]=$end if(!$h2{$g}->[1] || $end > $h2{$g}->[1]); $h2{$g}->[2]=$ch; $h2{$g}->[3]=$strand; $len=($end-$start)+1; $f=~/transcript_id "([^"]+)"/; $t=$1; push(@{$h{$g}->{$t}}, $len); END { for $g (keys %h) { ($st,$en,$chrm,$str)=@{$h2{$g}}; $sum=0; $min=(2**32)-1; $max=0; $mine=(2**32)-1; $maxe=0; $counte=0; $count=0; map { $i=0; for $e (@{$h{$g}->{$_}}) { $i+=$e; $mine=$e if($e < $mine); $maxe=$e if($e > $maxe); $counte++; } $count++; $sum+=$i; $max=$i if($i > $max); $min=$i if($i < $min); } (keys %{$h{$g}}); $st--; printf("$chrm\t$st\t$en\t$g\t$count\t$str\t$sum\t$min\t$max\t%.3f\t$counte\t$mine\t$maxe\t%.3f\n",($sum/$count),($sum/$counte));}}' | sort -k1,1 -k2,2n -k3,3n > ${ANNOTATION}.exons.stats.bed
