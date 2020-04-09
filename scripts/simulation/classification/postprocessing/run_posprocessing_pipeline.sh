#!/usr/bin/env bash
set -o pipefail -o nounset -o errexit 
#this runs the various scripts involved in post-processing
#predictions of the trained RF (both fl & nonfl) on real data
#it finds the predictions which agree with the trained class assignments
#and then further looks for intersections of those correct prediction regions
#across fl & nonfl sets to narrow down to the most confident predictions

BAM_BED_FILE=SRR1163655.sorted.bam.bed
CHROM_SIZES=SRR1163655.sorted.bam.chr_sizes

#e.g. SRR1163655.sorted.bam.bed.rl.nX3.minX2.mq.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm.sdX2.lmX4.lsX10
#or features.full as a symlink to the above (input to the prediction)
real_features=$1
fl_training_features=$2
nonfl_training_features=$3

#this script produces a "features.classes.bgz" file
/bin/bash -x join_features_and_classes.sh $real_features

/bin/bash -x extract_all_classes_from_training.sh $fl_training_features fl
/bin/bash -x extract_all_classes_from_training.sh $nonfl_training_features nonfl

#get predictions' recall & precision
#input prefixes to the 2 files, fl & nonfl, from the training extraction
/bin/bash -x pull_class_all_from_real_run.sh fl.all_classes.all.bed.sorted nonfl.all_classes.all.bed.sorted

#find agreement between fl & nonfl matched predictions (which agree with the training classes)
#use fl as the base set (more strict)
#also creates the BigBed (BB) files to use with the UCSC GB to visualize the output reads
/bin/bash -x find_fl_and_nonfl_agreement_matches.sh all.matches fl $BAM_BED_FILE $CHROM_SIZES

#finally run this at the root dir which has both "oxford" and "pacbio" as subdirs
#this will further find region intersections from the output set for each class from above
#but across technologies (pacbio & nonpore)
#currently skips problem-free since it's the largest category and takes the longest time
#even without it, this takes ~45min on stingray
/bin/bash -x join_ont_pb_categories.sh
