#!/bin/bash
set -o pipefail -o nounset -o errexit 
#creates list of exon coordinates + additional info useful in determining/ranking quality/novelty
BEDTOOLS='/data/bedtools2/bin/bedtools'
SWEEP='python ../sweep.py'
#params:

#1: file with original list of spliced exons extracted from BAM, format:
#chrom, start, end, strand, #_supporting_reads, comma-delimited list of supporting read indices, in_pacbio

#2: file with boundary coordinates of all long reads in original BAM (e.g. NA12878-DirectRNA.raw.cuff.transcript.coords), format:
#chrom, start, end, strand

#3: file with exon annotated coordinates in BED format 
#(e.g. recount_hg38_gencode_disjoint_exons.sorted.bed)

#4: file with gene annotated coordinates collapsed to be disjoint (no overlapping genes), format:
#chrom, start, end, etc...
#(e.g. gene_annotation_hg38/gencodev25_recount_genes.collapsed.sorted.tsv)

#5: file with list of pacbio exons (same format ONT exons in parameter 1)

#6: file with list of fully novel junctions from the same BAM (checked against annotation & Snaptron short reads)
#e.g. original_plus_4_wiggle.juncs.sorted.full.novel_1

#remove ones with less than 2 reads supporting
cat $1 | perl -ne 'chomp; $f=$_; ($c,$s,$e,$strand,$nr,$reads)=split(/\t/,$_); next if($nr < 2); print "$f\n";' | sort -k1,1 -k2,2n -k3,3n > ${1}.filtered
mv ${1} ${1}.old
ln -fs ${1}.filtered ${1}

#sort and determine which are matching pacbio exons
cat ${1} | perl -ne 'BEGIN { %h; open(IN,"<'${5}'"); while($line=<IN>) { chomp($line); ($c,$s,$e)=split(/\t/,$line); $k=join(":",($c,$s,$e)); $h{$k}=1; } close(IN); } chomp; $f=$_; ($c,$s,$e)=split(/\t/,$f); $k=join(":",($c,$s,$e)); $p=0; $p=1 if($h{$k}); print "$f\t$p\n";' > ${1}.pb

#collapse overlapping exons (required by later "sweep" step)
cat ${1}.pb | perl -ne 'chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads,$pac)=split(/\t/,$f); if($pc && $pc eq $c && $s <= $pe) { $p.=";".join("|",($s,$e,$o,$nr,$reads,$pac)); $pe = $e if($e > $pe); next; } if($pc) { print "$pc\t$ps\t$pe\t$p\n"; } $pc=$c; $ps=$s; $pe=$e; $p=join("|",($s,$e,$o,$nr,$reads,$pac)); END { if($pc) { print "$pc\t$ps\t$pe\t$p\n"; } }' > ${1}.pb.s.collasped

#get approximate count of *overlapping* long reads in the neighborhood (not exactly overlapping every exon due to collapsing)
$SWEEP --db-file <(sort -k1,1 -k2,2n -k3,3n ${2}) --query-file ${1}.pb.s.collasped --db-start-col 0 --q-start-col 0 --no-gzip --q-header --db-count > ${1}.pb.collapsed.lrsweep 

cat ${1}.pb.collapsed.lrsweep | perl -ne 'chomp; $f=$_; ($c_,$s_,$e_,$splices,$trcount)=split(/\t/,$f); @splices=split(/;/,$splices); for $splice (@splices) { @s=split(/\|/,$splice); ($s,$e)=splice(@s,0,2); $p=join("\t",@s); print "$c_\t$s\t$e\t$p\t$trcount\n";}' | sort -k1,1 -k2,2n -k3,3n > ${1}.pb.lroverlaps.s

#do this using BED files (start-1)
$BEDTOOLS closest -a <(cat ${1}.pb.lroverlaps.s | perl -ne 'chomp; @f=split(/\t/,$_); ($c,$s,$e)=splice(@f,0,3); $s--; print "$c\t$s\t$e\t".join("\t",@f)."\n";')  -b ${3} -d | perl -ne '$line=$_; chomp($line); @f=split(/\t/,$line); ($c,$s,$e,$o,$nr,$reads,$pac,$nlr)=splice(@f,0,8); $s++; $dist=pop(@f); $f[1]++; $f[0]=$f[3]; splice(@f,3,2); $rest=join("|",@f); $k=join(":",($c,$s,$e)); $h{$k}.=",$rest" if($h{$k}); $h{$k}="".join("\t",($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist))."\t".$rest if(!$h{$k}); END { for $v (values %h) { print "$v\n"; }}' | sort -k1,1 -k2,2n -k3,3n > ${1}.pb.lroverlaps.closest_merged.s

#further update "closest" distances to include annotated exons whose ends are far from the junction ends
#this allows for categorizing junctions which skip annotated exons and which result in clear alt. 5'/3' exon ends
cat ${1}.pb.lroverlaps.closest_merged.s | perl -ne 'BEGIN { $MAX=2**32; } chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$exons)=split(/\t/,$f); if($dist != 0) { print "$f\t".(($e-$s)+1)."\n"; next; } $min_dist = $MAX; for $ex (split(/,/,$exons)) { ($gid,$es,$ee,$eo)=split(/\|/,$ex); for $end (($es,$ee)) { for $end2 (($s,$e)) { $d=abs($end-$end2); $min_dist=$d if($d < $min_dist); }}} $dist = $min_dist; print "".join("\t",($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$exons,(($e-$s)+1)))."\n";' | egrep -v -e '(_alt)|(_random)|(_decoy)|(chrEBV)' > ${1}.pb.lroverlaps.closest_merged.s.inner_distance

#now filter for everything that's still within gene boundaries (but not overlapping an exon)
#$SWEEP --db-file ${1}.pb.lroverlaps.closest_merged.s --query-file ${4} --db-start-col 0 --q-start-col 0 --no-gzip | perl -ne 'chomp; @f=split(/\t/,$_); ($c,$s,$e,$info)=splice(@f,0,4); for $f (@f) { ($c2,$s2,$e2)=split(/\|\|\|/,$f); if($c eq $c2 && $s <= $s2 && $e >= $e2) { $f=~s/\|\|\|/\t/g; print "$f\n";}}' | sort -k1,1 -k2,2n -k3,3n > ${1}.pb.lroverlaps.closest_merged.overlapping_genes.s

#ln -fs ${1}.pb.lroverlaps.closest_merged.overlapping_genes.s pre_novel_exons
ln -fs ${1}.pb.lroverlaps.closest_merged.s.inner_distance pre_novel_exons


WIGGLE=11
cat pre_novel_exons | perl -ne 'BEGIN { $w='${WIGGLE}'; } chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$exons,$len)=split(/\t/,$f); $s1=$s-($w+1); $s2=$s+$w; $e1=$e-($w+1); $e2=$e+$w; print "".join("\t",($c,$s1,$s2,$o,$nr,$reads,$pac,$nlr,$dist,"$s-$e",0,$exons,$len))."\n"; print "".join("\t",($c,$e1,$e2,$o,$nr,$reads,$pac,$nlr,$dist,"$s-$e",1,$exons,$len))."\n";' | sort -k1,1 -k2,2n -k3,3n > pre_novel_exons.split_ends.${WIGGLE}.bed

#find non-overlaps with list of all short read exons looking for ends separately (to allow for distant containments/overlaps)
$BEDTOOLS intersect -sorted -a pre_novel_exons.split_ends.${WIGGLE}.bed -b <(zcat ../gtex_sra_junctions.split.bed.bgz) -v > novel_exons_bt.w${WIGGLE}.raw

cat novel_exons_bt.w${WIGGLE}.raw | perl -ne 'BEGIN { %h; } chomp; $f=$_; ($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$coord,$t,$exons,$len)=split(/\t/,$f); $tnot=!$t; $k=$coord."-".$tnot; if($h{$k}) { print "".$h{$k}."\n"; delete($h{$k}); next; } $k=$coord."-".$t; $coord=~/^(\d+)-(\d+)$/; $c1=$1; $c2=$2; $h{$k}=join("\t",($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$exons,$len));'  > novel_exons_bt.w${WIGGLE}.both_ends

cat novel_exons_bt.w${WIGGLE}.raw | perl -ne 'BEGIN { %h; } chomp; $f=$_; ($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$coord,$t,$exons,$len)=split(/\t/,$f); $k=$coord; $coord=~/^(\d+)-(\d+)$/; $c1=$1; $c2=$2; if(!$h{$k}) { print "".join("\t",($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$exons,$len))."\n"; } $h{$k}=1;' > novel_exons_bt.w${WIGGLE}.either_end

#finally, filter by novel junctions on both sides (exact):
cat ../gtex_junctions.tsv ../sra_junctions.tsv |  perl -ne 'BEGIN { %h; open(IN,"<pre_novel_exons"); $i=0; while($line=<IN>) { chomp($line); push(@a,$line); ($c,$s,$e)=split(/\t/,$line); $k1="$c:$s:1"; $k2="$c:$e:2"; $h{$k1}=1; $h{$k2}=1; } close(IN); } chomp; $f=$_; ($c,$s,$e)=split(/\t/,$f); $s--; $e++; $k1="$c:$s:2"; $k2="$c:$e:1"; delete($h{$k1}); delete($h{$k2}); END { for $a (@a) { ($c,$s,$e)=split(/\t/,$a); $k1="$c:$s:1"; $k2="$c:$e:2"; if($h{$k1} && $h{$k2}) { print "$a\n"; }}}' > novel_exons.exact.no_extra

cat novel_exons_bt.w10.either_end | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[8] == 0); if($f[6]==1) { print "$f\n"; } else { print STDERR "$f\n";}' > novel_exons_bt.w${WIGGLE}.either_end.pb 2> novel_exons_bt.w${WIGGLE}.either_end.npb

cat novel_exons.exact.no_extra | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[8] == 0); if($f[6]==1) { print "$f\n"; } else { print STDERR "$f\n";}' > novel_exons.exact.pb 2>novel_exons.exact.npb
