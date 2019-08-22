cut -f 2,3,5,6 gffcmp.*.gtf.tmap | egrep -e '	[ck=]	' | egrep -v -e '	1$' > tmap.multi_matching_contained_queries
wc -l tmap.multi_matching_contained_queries | wc -l
cut -f 3 tmap.multi_matching_contained_queries | sort -u > tmap.multi_matching_contained_queries.qids
wc -l tmap.multi_matching_contained_queries.qids

cut -f 2,3,4 gffcmp.*.gtf.refmap | egrep -e '	[ck=]	' > tmap.matching_contained_refs
fgrep -f tmap.multi_matching_contained_queries.qids.tabs tmap.matching_contained_refs > tmap.matching_contained_refs.multi_qids_filtered
wc -l tmap.multi_matching_contained_refs | wc -l
cut -f 3 tmap.multi_matching_contained_refs | sort -u > tmap.multi_matching_contained_refs.rids
wc -l tmap.multi_matching_contained_queries.rids
