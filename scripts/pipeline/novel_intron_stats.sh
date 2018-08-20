#!/bin/bash
set -o pipefail -o nounset -o errexit 
#creates list of junction coordinates + additional info useful in determining/ranking quality/novelty
BEDTOOLS='/data/bedtools2/bin/bedtools'
SWEEP='python ../../sweep.py'
#params:

#1: file with original list of spliced junctions extracted from BAM, format:
#chrom, start, end, strand, #_supporting_reads, comma-delimited list of supporting read indices

#2: file with boundary coordinates of all long reads in original BAM (e.g. NA12878-DirectRNA.raw.cuff.transcript.coords), format:
#chrom, start, end, strand

#3: file with junction annotated coordinates in BED format 
#(e.g. recount_hg38_gencode_disjoint_exons.sorted.bed)

#4: file with gene annotated coordinates collapsed to be disjoint (no overlapping genes), format:
#chrom, start, end, etc...
#(e.g. gene_annotation_hg38/gencodev25_recount_genes.collapsed.sorted.tsv)

#5: file with list of pacbio junctions (same format ONT junctions in parameter 1)
#so we can flag those which appear in both ONT and PB

#remove ones with less than 2 reads supporting
cat $1 | perl -ne 'chomp; $f=$_; ($c,$s,$e,$strand,$nr,$reads)=split(/\t/,$_); next if($nr < 2); print "$f\n";' | sort -k1,1 -k2,2n -k3,3n > ${1}.filtered
mv ${1} ${1}.old
ln -fs ${1}.filtered ${1}

#sort and determine which are matching pacbio junctions
cat ${1} | perl -ne 'BEGIN { %h; open(IN,"<'${5}'"); while($line=<IN>) { chomp($line); ($c,$s,$e)=split(/\t/,$line); $k=join(":",($c,$s,$e)); $h{$k}=1; } close(IN); } chomp; $f=$_; ($c,$s,$e)=split(/\t/,$f); $k=join(":",($c,$s,$e)); $p=0; $p=1 if($h{$k}); print "$f\t$p\n";' > ${1}.pb

#collapse overlapping junctions (required by later "sweep" step)
cat ${1}.pb | perl -ne 'chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads,$pac)=split(/\t/,$f); if($pc && $pc eq $c && $s <= $pe) { $p.=";".join("|",($s,$e,$o,$nr,$reads,$pac)); $pe = $e if($e > $pe); next; } if($pc) { print "$pc\t$ps\t$pe\t$p\n"; } $pc=$c; $ps=$s; $pe=$e; $p=join("|",($s,$e,$o,$nr,$reads,$pac)); END { if($pc) { print "$pc\t$ps\t$pe\t$p\n"; } }' > ${1}.pb.s.collasped

#get approximate count of *overlapping* long reads in the neighborhood (not exactly overlapping every junction due to collapsing)
$SWEEP --db-file <(sort -k1,1 -k2,2n -k3,3n ${2}) --query-file ${1}.pb.s.collasped --db-start-col 0 --q-start-col 0 --no-gzip --q-header --db-count > ${1}.pb.collapsed.lrsweep 
#$SWEEP --db-file <(sort -k1,1 -k2,2n -k3,3n ${2}) --query-file ${1}.pb.s.collasped --db-start-col 0 --q-start-col 0 --no-gzip --q-header > ${1}.pb.collapsed.lrsweep 

#cat ${1}.pb.collapsed.lrsweep | perl -ne 'chomp; $f=$_; @transcripts=split(/\t/,$f); ($c_,$s_,$e_,$splices)=splice(@transcripts,0,4); @splices=split(/;/,$splices); for $splice (@splices) { @s=split(/\|/,$splice); ($s,$e)=splice(@s,0,2); $p=join("\t",@s); $tcount=0; for $t (@transcripts) { ($tc,$ts,$te,$to)=split(/\|\|\|/,$t); if($c_ eq $tc && $s >= $ts && $e <= $te) { $tcount++; }} print "$c_\t$s\t$e\t$p\t$tcount\n";}' | sort -k1,1 -k2,2n -k3,3n > ${1}.pb.lroverlaps.s2
cat ${1}.pb.collapsed.lrsweep | perl -ne 'chomp; $f=$_; ($c_,$s_,$e_,$splices,$trcount)=split(/\t/,$f); @splices=split(/;/,$splices); for $splice (@splices) { @s=split(/\|/,$splice); ($s,$e)=splice(@s,0,2); $p=join("\t",@s); print "$c_\t$s\t$e\t$p\t$trcount\n";}' | sort -k1,1 -k2,2n -k3,3n > ${1}.pb.lroverlaps.s

#get distance to exons, do this using BED files (start-1), bedtools reports intervals in -a first
$BEDTOOLS closest -a <(cat ${1}.pb.lroverlaps.s | perl -ne 'chomp; @f=split(/\t/,$_); ($c,$s,$e)=splice(@f,0,3); $s--; print "$c\t$s\t$e\t".join("\t",@f)."\n";')  -b ${3} -d | perl -ne '$line=$_; chomp($line); @f=split(/\t/,$line); ($c,$s,$e,$o,$nr,$reads,$pac,$nlr)=splice(@f,0,8); $s++; $dist=pop(@f); $f[1]++; $f[0]=$f[3]; splice(@f,3,2); $rest=join("|",@f); $k=join(":",($c,$s,$e)); $h{$k}.=",$rest" if($h{$k}); $h{$k}="".join("\t",($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist))."\t".$rest if(!$h{$k}); END { for $v (values %h) { print "$v\n"; }}' | sort -k1,1 -k2,2n -k3,3n > ${1}.pb.lroverlaps.closest_merged.s

#further update "closest" distances to include annotated exons whose ends are far from the junction ends
#this allows for categorizing junctions which skip annotated exons and which result in clear alt. 5'/3' exon ends
cat ${1}.pb.lroverlaps.closest_merged.s | perl -ne 'BEGIN { $MAX=2**32; } chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$exons)=split(/\t/,$f); if($dist != 0) { print "$f\n"; next; } $min_dist = $MAX; for $ex (split(/,/,$exons)) { ($gid,$es,$ee,$eo)=split(/\|/,$ex); for $end (($es,$ee)) { for $end2 (($s,$e)) { $d=abs($end-$end2); $min_dist=$d if($d < $min_dist); }}} $dist = $min_dist; print "".join("\t",($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$exons))."\n";' > ${1}.pb.lroverlaps.closest_merged.s.inner_distance

#now filter for everything that's still within gene boundaries (but not overlapping an exon)
#UPDATE: don't limit to those within gene boundaries
#$SWEEP --db-file ${1}.pb.lroverlaps.closest_merged.s --query-file ${4} --db-start-col 0 --q-start-col 0 --no-gzip | perl -ne 'chomp; @f=split(/\t/,$_); ($c,$s,$e,$info)=splice(@f,0,4); for $f (@f) { ($c2,$s2,$e2)=split(/\|\|\|/,$f); if($c eq $c2 && $s <= $s2 && $e >= $e2) { $f=~s/\|\|\|/\t/g; print "$f\n";}}' | sort -k1,1 -k2,2n -k3,3n > ${1}.pb.lroverlaps.closest_merged.overlapping_genes.s

#ln -fs ${1}.pb.lroverlaps.closest_merged.overlapping_genes.s pre_novel_junctions
ln -fs ${1}.pb.lroverlaps.closest_merged.s.inner_distance pre_novel_junctions

#sort and collapse overlapping junctions (required by later "sweep" step), also bump the start/ends by -/+ 10 to account for wiggle
#cat pre_novel_junctions | perl -ne 'chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$junctions)=split(/\t/,$f); $s-=10; $e+=10; if($pc && $pc eq $c && $s <= $pe) { $p.=";".join("||",($s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$junctions)); $pe = $e if($e > $pe); next; } if($pc) { print "$pc\t$ps\t$pe\t$p\n"; } $pc=$c; $ps=$s; $pe=$e; $p=join("||",($s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$junctions)); END { if($pc) { print "$pc\t$ps\t$pe\t$p\n"; } }' > pre_novel_junctions.collapsed

#just take coords and strand and then track by either (both) ends
WIGGLE=10
cat pre_novel_junctions | perl -ne 'BEGIN { $w='${WIGGLE}'; } chomp; $f=$_; ($c,$s,$e,$o,$nr,$reads,$pac,$nlr,$dist,$junctions)=split(/\t/,$f); $s1=$s-($w+1); $s2=$s+$w; $e1=$e-($w+1); $e2=$e+$w; print "".join("\t",($c,$s1,$s2,$o,$nr,$reads,$pac,$nlr,$dist,"$s-$e",0,$junctions))."\n"; print "".join("\t",($c,$e1,$e2,$o,$nr,$reads,$pac,$nlr,$dist,"$s-$e",1,$junctions))."\n";' | sort -k1,1 -k2,2n -k3,3n > pre_novel_junctions.split_ends.${WIGGLE}.bed

#find non-overlaps with list of all short read junctions looking for ends separately (to allow for distant containments/overlaps)
$BEDTOOLS intersect -sorted -a pre_novel_junctions.split_ends.${WIGGLE}.bed -b <(zcat ../../gtex_sra_junctions.split.bed.bgz) -v | perl -ne 'BEGIN { %h; } chomp; $f=$_; ($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$cord,$t,$junctions)=split(/\t/,$f); $tnot=!$t; $k=$coord."-".$tnot; if($h{$k}) { print "".$h{$k}."\n"; delete($h{$k}); next; } $k=$coord."-".$t; $coord=~/^(\d+)-(\d+)$/; $c1=$1; $c2=$2; $h{$k}=join("\t",($c,$c1,$c2,$o,$nr,$reads,$pac,$nlr,$dist,$junctions));' > novel_junctions_bedtools.wiggle.${WIGGLE}

#finally filter to get novel junctions w/ wiggle
#$SWEEP --db-file <(zcat ../../gtex_sra_junctions.tsv.gz) --query-file pre_novel_junctions.collapsed --db-start-col 0 --q-start-col 0 --no-gzip --not-matched > novel_junctions.wiggle

#filter out uninteresting contigs
egrep -v -e '(_alt)|(_random)|(_decoy)' novel_junctions_bedtools.wiggle.${WIGGLE} > novel_junctions_bedtools.wiggle.${WIGGLE}.no_alts_randoms_decoys

#finally, filter by novel junctions w/o wiggle (exact):
cat ../../gtex_junctions.tsv ../../sra_junctions.tsv |  perl -ne 'BEGIN { %h; open(IN,"<pre_novel_junctions"); $i=0; while($line=<IN>) { chomp($line); push(@a,$line); ($c,$s,$e)=split(/\t/,$line); $k="$c:$s:$e"; $h{$k}=1; } close(IN); } chomp; $f=$_; ($c,$s,$e)=split(/\t/,$f); $k="$c:$s:$e"; delete($h{$k}); END { for $a (@a) { ($c,$s,$e)=split(/\t/,$a); $k="$c:$s:$e"; if($h{$k}) { print "$a\n"; }}}' > novel_junctions.exact

egrep -v -e '(_alt)|(_random)|(_decoy)' novel_junctions.exact > novel_junctions.exact.no_alts_randoms_decoys

#ln -fs ${1}.pb.lroverlaps.closest_merged.overlapping_genes.s.novel novel_junctions
#sort -k9,9nr -k7,7nr -k5,5nr novel_junctions > novel_junctions.sorted
cat novel_junctions_bedtools.wiggle.${WIGGLE}.no_alts_randoms_decoys | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[8] == 0); if($f[6]==1) { print "$f\n"; } else { print STDERR "$f\n";}' > novel_junctions_bedtools.wiggle.${WIGGLE}.pb 2> novel_junctions_bedtools.wiggle.${WIGGLE}.npb

cat novel_junctions.exact.no_alts_randoms_decoys | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); next if($f[8] == 0); if($f[6]==1) { print "$f\n"; } else { print STDERR "$f\n";}' > novel_junctions.exact.pb 2>novel_junctions.exact.npb
