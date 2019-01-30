#!/usr/bin/env bash
#takes the output of find_fl_and_nonfl_agreement_matches.sh
#and joins between oxford and pacbio based on coordinate overlap
#assumes we're run from the parent directory where
#"oxford" and "pacbio" are subdirs

o="./oxford"
p="./pacbio"
suffix="all.matches.joined"

#for f in problem-free non-recurrent recurrent novel; do
for f in non-recurrent recurrent novel; do
    o1="${o}/${f}.${suffix}"
    cat $o1 | perl -ne 'chomp; $f=$_; $f=~s/:/\t/; $f=~s/-/\t/; print "$f\n";' | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n > ${o1}.tab.sorted
    bedtools intersect -sorted -wo -a ${o1}.tab.sorted -b ${o1}.tab.sorted | perl -ne 'chomp; ($c,$s,$e,$c2,$s2,$e2)=split(/\t/,$_); $e1="$c:$s-$e"; $e2="$c2:$s2-$e2"; next if($h{$e1}); $h{$e2}=1 if($e1 ne $e2); print "$e1\n";' | sort -u > ${f}.ont.intersecting.coords
    p1="${p}/${f}.${suffix}"
    cat $p1 | perl -ne 'chomp; $f=$_; $f=~s/:/\t/; $f=~s/-/\t/; print "$f\n";' | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n > ${p1}.tab.sorted
    bedtools intersect -sorted -wo -a ${p1}.tab.sorted -b ${p1}.tab.sorted | perl -ne 'chomp; ($c,$s,$e,$c2,$s2,$e2)=split(/\t/,$_); $e1="$c:$s-$e"; $e2="$c2:$s2-$e2"; next if($h{$e1}); $h{$e2}=1 if($e1 ne $e2); print "$e1\n";' | sort -u > ${f}.pb.intersecting.coords
    bedtools intersect -sorted -wo -a <(cat ${f}.ont.intersecting.coords| perl -ne 'chomp; $f=$_; $f=~s/:/\t/; $f=~s/-/\t/; print "$f\n";' | sort -k1,1 -k2,2n -k3,3n) -b <(cat ${f}.pb.intersecting.coords | perl -ne 'chomp; $f=$_; $f=~s/:/\t/; $f=~s/-/\t/; print "$f\n";' | sort -k1,1 -k2,2n -k3,3n) | perl -ne 'chomp; ($c,$s,$e,$c2,$s2,$e2)=split(/\t/,$_); $e1="$c:$s-$e"; $e2="$c2:$s2-$e2"; next if($h{$e1}); $h{$e2}=1 if($e1 ne $e2); print "$e1\n";' | sort -u > ${f}.ont.pb.intersecting.coords
done
