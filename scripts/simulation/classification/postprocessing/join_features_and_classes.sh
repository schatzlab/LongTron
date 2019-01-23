#!/usr/bin/env bash
features_f=$1
p='features.full.just_features'

#paste $features_f ${p}.nonfl.2_class_out ${p}.fl.2_class_out ${p}.nonfl.4_class_out ${p}.fl.4_class_out | bgzip > features.classes.bgz
paste $features_f ${p}.nonfl.5_class_out ${p}.fl.5_class_out | bgzip > features.classes.bgz
tabix -s1 -b2 -e3 features.classes.bgz

