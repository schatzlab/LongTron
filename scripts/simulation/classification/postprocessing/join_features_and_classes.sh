#!/usr/bin/env bash
#this pastes together the input features from the real data
#with the RF class probabilities for both fl & nonfl ("5 class")
#only use the probabilities from here on

#e.g. SRR1163655.sorted.bam.bed.rl.nX3.minX2.mq.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm.sdX2.lmX4.lsX10
#or features.full as a symlink to the above (input to the prediction)
features_f=$1
p='features.full.just_features'

#paste $features_f ${p}.nonfl.2_class_out ${p}.fl.2_class_out ${p}.nonfl.4_class_out ${p}.fl.4_class_out | bgzip > features.classes.bgz
paste $features_f ${p}.nonfl.5_class_out ${p}.fl.5_class_out | bgzip > features.classes.bgz
tabix -s1 -b2 -e3 features.classes.bgz

