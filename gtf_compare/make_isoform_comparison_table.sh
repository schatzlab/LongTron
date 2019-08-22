
#/path/to/<comparison>.tmap file
tmap=$1

#cut -f 5 $tmap | sort -u > ${tmap}.unique_qs
#total_q=`wc -l ${tmap}.unique_qs | cut -d' ' -f 1`

cut -f 3 $tmap | sort | uniq -c > ${tmap}.counts
#only want 1) multi-exon 2) full matching "=", contains matching ("c" for contained *in* ref, and "k" for contained in query)
cut -f 2,3,5,6 $tmap | fgrep "	=	" | egrep -v -e '	1$' > ${tmap}.exact.multi
#query is contained in ref
cut -f 2,3,5,6 $tmap | fgrep "	c	" | egrep -v -e '	1$' > ${tmap}.contained.multi
#query contains ref
cut -f 2,3,5,6 $tmap | fgrep "	k	" | egrep -v -e '	1$' > ${tmap}.contains.multi

#cat ${tmap}.exact.multi ${tmap}.contains.multi ${tmap}.contained.multi > ${tmap}.both.multi

#now get intergenic, completely novel multi-exon queries
cut -f 2,3,5,6 $tmap | fgrep "	u	" | egrep -v -e '	1$' > ${tmap}.novel.multi
cut -f 2,3,5,6 $tmap | fgrep "	p	" | egrep -v -e '	1$' >> ${tmap}.novel.multi

#cut -f 5 ${tmap}.novel.multi | sort -u > ${tmap}.novel.multi.unique_qs

#repeat regions
cut -f 2,3,5,6 $tmap | fgrep "	r	" | egrep -v -e '	1$' > ${tmap}.repeats.multi
#cut -f 5 ${tmap}.repeats.multi | sort -u > ${tmap}.repeats.multi.unique_qs

#cut -f 5 ${tmap}.both.multi | sort -u > ${tmap}.both.multi.unique_qs
total_q=`cut -f 5 $tmap | sort -u | wc -l | cut -d' ' -f 1`
matching_q=`cut -f 3 ${tmap}.exact.multi | sort -u | wc -l | cut -d' ' -f 1`
contained_q=`cut -f 3 ${tmap}.contained.multi | sort -u | wc -l | cut -d' ' -f 1`
contains_q=`cut -f 3 ${tmap}.contains.multi | sort -u | wc -l | cut -d' ' -f 1`
novel_q=`cut -f 3 ${tmap}.novel.multi | sort -u | wc -l | cut -d' ' -f 1`
repeats_q=`cut -f 3 ${tmap}.repeats.multi | sort -u | wc -l | cut -d' ' -f 1`

#print out: <dataset>	all_matching_contained_contains,all_matching,contained,contains,non_matching_overlaps,novel,repeats
dataset=`echo "$tmap" | perl -ne 'chomp; $f=$_; $t='$total_q'; $m='$matching_q'; $cd='$contained_q'; $cs='$contains_q'; $n='$novel_q'; $r='$repeats_q'; @f=split(/\//,$f); $f=pop(@f); $dataset=pop(@f); $mp=100*($m/$t); $cdp=100*($cd/$t); $csp=100*($cs/$t); $mcds=$m+$cd+$cs; $mcdsp=100*($mcds/$t); $np=100*($n/$t); $rp=100*($r/$t); $o=$t-($m+$cd+$cs+$n+$r); $op=100*($o/$t); printf("$dataset\t$t\t%.1f\%,%.1f\%,%.1f\%,%.1f\%,%.1f\%,%.1f\%,%1.f\%\t$mcds,$m,$cd,$cs,$o,$n,$r\n",$mcdsp,$mp,$cdp,$csp,$op,$np,$rp);'`
echo $dataset
