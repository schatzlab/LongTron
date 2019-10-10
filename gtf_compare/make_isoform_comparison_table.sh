
#/path/to/<comparison>.tmap file
tmap=$1

cut -f 2,3,5,6 $tmap | egrep -v -e '	1$' > ${tmap}.multi
cut -f 2 ${tmap}.multi | sort | uniq -c > ${tmap}.counts


egrep -v -e '	(=|c|k|u|p|r)	' ${tmap}.multi > ${tmap}.wrong.multi
cut -f 2 ${tmap}.wrong.multi | sort | uniq -c > ${tmap}.wrong.multi.counts

#only want 1) multi-exon 2) full matching "=", contains matching ("c" for contained *in* ref, and "k" for contained in query)
fgrep "	=	" ${tmap}.multi > ${tmap}.exact.multi
#query is contained in ref
fgrep "	c	" ${tmap}.multi > ${tmap}.contained.multi
#query contains ref
fgrep "	k	" ${tmap}.multi > ${tmap}.contains.multi

#now get intergenic, completely novel multi-exon queries
fgrep "	u	" ${tmap}.multi > ${tmap}.novel.multi
fgrep "	p	" ${tmap}.multi >> ${tmap}.novel.multi

#repeat regions
fgrep "	r	" ${tmap}.multi > ${tmap}.repeats.multi

#get total count of all *multi-exon* query isofrags
total_q=`cut -f 3 ${tmap}.multi | tail -n+2 | sort -u | wc -l | cut -d' ' -f 1`
matching_q=`cut -f 3 ${tmap}.exact.multi | sort -u | wc -l | cut -d' ' -f 1`
contained_q=`cut -f 3 ${tmap}.contained.multi | sort -u | wc -l | cut -d' ' -f 1`
contains_q=`cut -f 3 ${tmap}.contains.multi | sort -u | wc -l | cut -d' ' -f 1`
novel_q=`cut -f 3 ${tmap}.novel.multi | sort -u | wc -l | cut -d' ' -f 1`
repeats_q=`cut -f 3 ${tmap}.repeats.multi | sort -u | wc -l | cut -d' ' -f 1`

#print out: <dataset>	all_matching_contained_contains,all_matching,contained,contains,non_matching_overlaps,novel,repeats
dataset=`echo "$tmap" | perl -ne 'chomp; $f=$_; $t='$total_q'; $m='$matching_q'; $cd='$contained_q'; $cs='$contains_q'; $n='$novel_q'; $r='$repeats_q'; @f=split(/\//,$f); $f=pop(@f); $dataset=pop(@f); $mp=100*($m/$t); $cdp=100*($cd/$t); $csp=100*($cs/$t); $mcds=$m+$cd+$cs; $mcdsp=100*($mcds/$t); $np=100*($n/$t); $rp=100*($r/$t); $o=$t-($m+$cd+$cs+$n+$r); $op=100*($o/$t); printf("$dataset\t$t\t%.1f\%,%.1f\%,%.1f\%,%.1f\%,%.1f\%,%.1f\%,%1.f\%\t$mcds,$m,$cd,$cs,$o,$n,$r\n",$mcdsp,$mp,$cdp,$csp,$op,$np,$rp);'`
echo $dataset

