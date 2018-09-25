#!/bin/bash
#used to add sequence overlap features to a set of intervals (transcript, reads, etc...)
#this is primarily used for classification in terms of long read alignment error profiling

#Features:
#1) Repeat Masker overlap
#2) Tandem Repeat Finder (TRF) overlap (subset of #1)
#3) Splice motif count
#4) # of exons
#5) size of exons
#6) size of introns
#7) # of overlapping reads/transcripts
#8) # of common 150 SNPs overlapping

SNP=snp150Common.txt.gz
RM=hg38_repeatmasker_rmsk.gz

#assume we get a true BED file as input
IN=$1

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
