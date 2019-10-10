
ln -fs na12878_ox_vs_all na12878_oxford_vs_all
ln -fs na12878_pb_vs_all na12878_pacbio_vs_all

for f in oxford pacbio ; do
    cut -f 3 na12878_${f}_vs_all/gffcmp.${f}_na12878.gtf.tmap.novel.multi | sort -u | perl -ne 'chomp; print "transcript_id \"$_\"\n";' > ${f}.novel.tids
    fgrep -f ${f}.novel.tids ${f}_na12878.gtf > ${f}_na12878.novels.gtf
done

#ln -fs ../oxford_na12878.novels.gtf na12878_oxn_vs_pbn/oxford_na12878.novels.gtf
#ln -fs ../pacbio_na12878.novels.gtf na12878_oxn_vs_pbn/pacbio_na12878.novels.gtf
/bin/bash -x ./run_gcom.sh na12878_oxn_vs_pbn pacbio_na12878.novels.gtf oxford_na12878.novels.gtf 

#ln -fs ../pacbio_na12878.novels.gtf na12878_pbn_vs_oxn/pacbio_na12878.novels.gtf
#ln -fs ../oxford_na12878.novels.gtf na12878_pbn_vs_oxn/oxford_na12878.novels.gtf
/bin/bash -x ./run_gcom.sh na12878_pbn_vs_oxn oxford_na12878.novels.gtf pacbio_na12878.novels.gtf
