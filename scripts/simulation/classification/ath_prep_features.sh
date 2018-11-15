#!/bin/bash
#for non-human, where there are less resources
ANNOT_VER=TAIR10
ANNOTATION=Araport11_GFF3_genes_transposons.201606.gff
CHROM_SIZES=TAIR10_Chr.all.fasta.chrm_sizes
SNPS=TAIR9_GFF3_polymorphisms.gff
NUM_CHRMS=7
GENOME=TAIR10_Chr.all.fasta

cat ${GENOME} | perl -ne 'chomp; $f=$_; if($f=~/^>/) { ($j,$c)=split(/[>\s+]/,$f); if($pc) { print "$pc\t".(length($seq))."\n"; } $pc=$c; $seq=""; next; } $seq.=$f; END {if($pc) { print "$pc\t".(length($seq))."\n"; }}' > $CHROM_SIZES

#first fix Ath GFF to match what GFFREAD pulled out into FASTA
cat ${ANNOTATION} | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($type,$i)=($f[2],$f[8]); $i=~/ID=([^;]+);/; $id=$1; $h{$id}=$type; $i=~/Parent=([^;]+)/; $p=$1;  next if($type eq "transcript_region" || ($type eq "exon" && $h{$p} eq "transcript_region")); if($type=~/^miRNA(_primary_transcript)?/) { $f[2]=$type."_exon"; } print "".join("\t",@f)."\n";' > ${ANNOTATION}_

cat ${ANNOTATION}_ | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); if($f[2] eq "miRNA_exon") { $f=~s/\tID=([^;]+)/\tID=$1;Parent=$1/; } if($f[2] eq "miRNA_primary_transcript_exon") { $f=~s/Parent=/Derives_from=/; $f=~s/\tID=([^;]+)/\tID=$1;Parent=$1/; } print "$f\n";' > ${ANNOTATION}.matches_fasta

ANNOTATION=Araport11_GFF3_genes_transposons.201606.gff.matches_fasta

#pull out annotated junctions
egrep -e 'exon	'  $ANNOTATION | sort -k1,1 -k4,4n -k5,5n | cut -f 1,4,5,7,9 | perl -ne 'BEGIN { open(OUT,">'${ANNOTATION}'.exons.bed"); } chomp; $f=$_; ($c,$s,$e,$o,$info)=split(/\t/,$f); $info=~/Parent=([^;]+);/; $t=$1; @t=split(/,/,$t); $info=~/exon:(\d+)/; $eid=$1; for $t (@t) { $t.= "_0_$o"; print OUT "$c\t".($s-1)."\t$e\t$eid;$t\t0\t$o\n"; if($ts{$t}) { $pe=$ts{$t}; $c1=$pe+1; $c2=$s-1;  $eid--; print "$c\t$c1\t$c2\t$t\t$eid\n"; } $ts{$t}=$e;} END { close(OUT); }' | sort -k1,1 -k2,2n -k3,3n | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($c,$s,$e,$n)=@f; if($n=~/_\+$/) { print "$f\n"; $f=join("\t",@f); next; } if(!$ns{$n}) { $eid=0; } else { $eid=$ns{$n} } $eid++; print "$c\t$s\t$e\t$n\t$eid\n"; $ns{$n}=$eid;'  > ${ANNOTATION}.junctions

#get total set of transcripts from exons
cat ${ANNOTATION}.exons.bed | perl -ne 'chomp; ($c,$s,$e,$n,$junk,$o)=split(/\t/,$_); ($eid,$tid)=split(/;/,$n); if(!$t{$tid}) { $t{$tid}=[$c,$s,$e,$o]; } if($t{$tid}->[2] < $e) { $t{$tid}->[2]=$e; } END { for $tid (keys %t) { ($tc,$ts,$te,$to)=@{$t{$tid}}; print "$tc\t$ts\t$te\t$tid\t0\t$to\n"; }}' | sort -k1,1 -k2,2n -k3,3n > ${ANNOTATION}.transcripts.bed

#get number of exon per transcript
egrep -e 'exon	' $ANNOTATION | perl -ne 'chomp; $f=$_; $f=~/Parent=([^;]+);/; $t=$1; @t=split(/,/,$t); for $t (@t) { print "$t\n"; }' | sort | uniq -c | sort -k1,1nr | perl -ne 'chomp; ($j,$c,$n)=split(/\s+/,$_); print "$n\t$c\n";' > ${ANNOTATION}.transcripts_exon_count

#convert SNPs to BED
cut -f 1,4,5,9 $SNPS | perl -ne 'chomp; ($c,$s,$e,$info)=split(/\t/,$_); $info=~/ID=([^;]+);/; $n=$1; $s--; print "$c\t$s\t$e\t$n\n";' | sort -k1,1 -k2,2n -k3,3n > ${SNPS}.sorted.bed

#prep GC/UMAP track data by removing non-cannonical chromosomes and sort
#TODO: pull out GC content and mappability for Ath

#create perbase for exon/transcript densities
#this assumes that there are multiple overlapping exons in some regions (maybe even duplicates across transcripts)
#these are big files, so have the ANNOTATION file sitting on a large filesystem with plenty of free space
cat ${ANNOTATION}.exons.bed | sort -k1,1 -k2,2n -k3,3n | perl -ne 'BEGIN { open(OUT,">'${ANNOTATION}'.exons.sorted.bed"); } chomp; $f=$_; ($c,$s,$e,$name,$score,$o)=split(/\t/,$f); print OUT "$c\t$s\t$e\n"; if(!$pc || ($pc ne $c || $s > $pe)) { if($pc) { for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; }} $ps=$s; $pc=$c; $pe=$e; $po=$o; next; } $pe=$e; END { if($pc) { for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; } } close(OUT);}' | bgzip > ${ANNOTATION}.exons.perbase.bgz

bedtools intersect -sorted -c -a <(zcat ${ANNOTATION}.exons.perbase.bgz) -b <(sort -k1,1 -k2,2n -k3,3n ${ANNOTATION}.exons.sorted.bed) | bgzip > ${ANNOTATION}.exons.perbase.counts.bgz

cat ${ANNOTATION}.transcripts.bed | sort -k1,1 -k2,2n -k3,3n | perl -ne 'BEGIN { open(OUT,">'${ANNOTATION}'.transcripts.sorted.bed"); } chomp; $f=$_; ($c,$s,$e,$name,$score,$o)=split(/\t/,$f); print OUT "$c\t$s\t$e\n"; if(!$pc || ($pc ne $c || $s > $pe)) { if($pc) { for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; }} $ps=$s; $pc=$c; $pe=$e; $po=$o; next; } $pe=$e; END { if($pc) { for $i ($ps..($pe-1)) { print "$pc\t$i\t".($i+1)."\n"; } } close(OUT);}' | bgzip > ${ANNOTATION}.transcripts.perbase.bgz

bedtools intersect -sorted -c -a <(zcat ${ANNOTATION}.transcripts.perbase.bgz) -b <(sort -k1,1 -k2,2n -k3,3n ${ANNOTATION}.transcripts.sorted.bed) | bgzip > ${ANNOTATION}.transcripts.perbase.counts.bgz

#extract out all cannonical splice motifs (forward + reverese same for all eukaryotes) from reference across all main chromosomes
#takes ~33m22.814s for HG38
#takes ~2m for Ath
(head -${NUM_CHRMS} $CHROM_SIZES | cut -f 1 | sort | perl -ne 'BEGIN { %F=("GT"=>1,"gt"=>1,"Gt"=>1,"gT"=>1,"AG"=>1,"ag"=>1,"Ag"=>1,"aG"=>1); %R=("CT"=>1,"AC"=>1,"cT"=>1,"Ct"=>1,"ct"=>1,"Ac"=>1,"aC"=>1,"ac"=>1); } chomp; $chrm=$_; @s=`samtools faidx '${GENOME}' $chrm | fgrep -v ">"`; chomp(@s); $seq=join("",@s); $len=length($seq); for($i=0;$i<($len-1);$i++) { $m=substr($seq,$i,2); $b="$chrm\t$i\t".($i+1)."\t$m\n"; if($F{$m}) { print $b; } if($R{$m}) { print STDERR $b; } }' | bgzip > ${ANNOT_VER}.forward_splice_motifs.all.tsv.bgz) 2>&1 | bgzip > ${ANNOT_VER}.reverse_splice_motifs.all.tsv.bgz

#merge forward and reverse strands of all cannonical splice motifs from reference into one file (sorted)
cat <(zcat ${ANNOT_VER}.forward_splice_motifs.all.tsv.bgz | perl -ne 'chomp; print "$_\t0\t+\n";') <(zcat ${ANNOT_VER}.reverse_splice_motifs.all.tsv.bgz | perl -ne 'chomp; print "$_\t0\t-\n";') | sort -k1,1 -k2,2n -k3,3n | bgzip > ${ANNOT_VER}.splice_motifs.all.tsv.bgz


#get per-gene level stats about transcript and exon lengths/counts/sums/averages
cat  ${ANNOTATION} | egrep -e 'exon	' | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); ($ch,$start,$end,$strand)=($f[0],$f[3],$f[4],$f[6]); $f=~/\tID=([^:;]+)/; $g=$1; $h2{$g}->[0]=$start if(!$h2{$g}->[0] || $start < $h2{$g}->[0]); $h2{$g}->[1]=$end if(!$h2{$g}->[1] || $end > $h2{$g}->[1]); $h2{$g}->[2]=$ch; $h2{$g}->[3]=$strand; $len=($end-$start)+1; $f=~/Parent=([^;]+);/; $t=$1; push(@{$h{$g}->{$t}}, $len); END { for $g (keys %h) { ($st,$en,$chrm,$str)=@{$h2{$g}}; $sum=0; $min=(2**32)-1; $max=0; $mine=(2**32)-1; $maxe=0; $counte=0; $count=0; map { $i=0; for $e (@{$h{$g}->{$_}}) { $i+=$e; $mine=$e if($e < $mine); $maxe=$e if($e > $maxe); $counte++; } $count++; $sum+=$i; $max=$i if($i > $max); $min=$i if($i < $min); } (keys %{$h{$g}}); printf("$chrm\t$st\t$en\t$g\t$count\t$str\t$sum\t$min\t$max\t%.3f\t$counte\t$mine\t$maxe\t%.3f\n",($sum/$count),($sum/$counte));}}' | sort -k1,1 -k2,2n -k3,3n > ${ANNOTATION}.exons.stats.bed
