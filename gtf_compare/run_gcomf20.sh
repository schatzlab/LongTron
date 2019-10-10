#!/usr/bin/env bash
set -ex -o pipefail

dir=$1
gtf1=$2
gtf2=$3

mkdir -p $dir
cd $dir
ln -fs ../${gtf1}
ln -fs ../${gtf2}
../gcom -r $gtf1 --fuzz-length 20 $gtf2 > run 2>&1
/bin/bash -x ../../make_isoform_comparison_table.sh gffcmp.${gtf2}.tmap > comparison.tsv
#make backup of 1st run's generically named files
#if pacbio NA12878 flip the strands of those transcripts which ended up in O/S categories
FLIP_STRAND=`echo $dir | perl -ne '$n=$_; chomp($n); if($n=~/pb/) { print "1"; } else { print "0"; }'`
if [[ "$FLIP_STRAND" == "1" ]]; then
    egrep -e '	(o|s)	' gffcmp.${gtf2}.tmap.multi | cut -f 3 | sort -u  | sed -e 's/^\(.*\)$/"\1"/' > gffcmp.${gtf2}.tmap.multi.o_s
    cat gffcmp.${gtf2}.tmap.multi.o_s $gtf2 | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); if(scalar @f == 1) { $f=~s/"//g; $h{$f}=1; next; } $f[8]=~/transcript_id "([^"]+)"/; $tid=$1; if($h{$tid}) { $s=$f[6]; $f[6]=($s eq "+"?"-":"+"); } print "".join("\t",@f)."\n";' > ${gtf2}.o_s_strand_flipped
    mkdir ../${dir}.flipped
    cd ../${dir}.flipped
    ln -fs ../${gtf1}
    ln -fs ../${dir}/${gtf2}.o_s_strand_flipped
    #now re-run with flipped strands
    ../gcom -r $gtf1 --fuzz-length 20 ${gtf2}.o_s_strand_flipped > run 2>&1
    /bin/bash -x ../../make_isoform_comparison_table.sh gffcmp.${gtf2}.o_s_strand_flipped.tmap > comparison.tsv
fi
cd ../
