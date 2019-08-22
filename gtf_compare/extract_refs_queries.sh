fgrep "icF" run > run.icf
cut -f 2 run.icf | sort -u > run.icf.refs
cut -f 3 run.icf | sort -u > run.icf.qs

cut -f 3 *.tmap | sort | uniq -c > tmap.counts

cut -f 2,3,5,6 *.tmap | fgrep "	=	" | egrep -v -e '	1$' > tmap.exact.multi
cut -f 2,3,5,6 *.tmap | fgrep "	k	" | egrep -v -e '	1$' > tmap.contains.multi
cut -f 2,3,5,6 *.tmap | fgrep "	c	" | egrep -v -e '	1$' >> tmap.contains.multi

cat tmap.exact.multi tmap.contains.multi > tmap.both.multi
cut -f 1 tmap.both.multi | sort -u > tmap.both.multi.refs
cut -f 3 tmap.both.multi | sort -u > tmap.both.multi.qs

for f in refs qs; do
    comm -1 -2 tmap.both.multi.${f} run.icf.${f} > ${f}.shared.both
    comm -1 -3 tmap.both.multi.${f} run.icf.${f} > ${f}.only_in.run
    comm -2 -3 tmap.both.multi.${f} run.icf.${f} > ${f}.only_in.t_map
done

fgrep -f qs.only_in.t_map tmap.both.multi | cut -f 1 | sort -u | perl -ne 'chomp; print "$_\t\n";' > qs.only_in.t_map.refs
cat run.icf.refs | perl -ne 'chomp; print "$_\t\n";' > run.icf.refs.tabs
fgrep -f qs.only_in.t_map.refs run.icf.refs.tabs > qs.only_in.t_map.refs.in_run
fgrep -v -f qs.only_in.t_map.refs.in_run qs.only_in.t_map.refs > qs.only_in.t_map.refs.not_in_run
fgrep -f qs.only_in.t_map.refs.not_in_run tmap.both.multi > tmap.both.multi.missing_refs
cat qs.only_in.t_map | perl -ne 'chomp; print "$_\t\n";' > qs.only_in.t_map.tabs
fgrep -f qs.only_in.t_map.tabs tmap.both.multi.missing_refs > tmap.both.multi.missing_refs.missing_qs

cat tmap.both.multi.missing_refs.missing_qs | perl -ne 'chomp; ($r,$m,$q,$n)=split(/\t/,$_); `fgrep $r na12878_ox.flair.collapse.isoforms.gtf  | cut -f 1,4,5,7 > $r`; `cat $r > $r.vs.$q`; `wc -l $r >> $r.vs.$q`; `fgrep $q na12878_pb.flair.collapse.isoforms.gtf  | cut -f 1,4,5,7 > $q`; `cat $q >> $r.vs.$q`; `wc -l $q >> $r.vs.$q`; `diff $r $q >> $r.vs.$q`;'


mkdir non_run_qs_check
mv ENST00000* non_run_qs_check/
mv SRR1163655.* non_run_qs_check/
mv *:* non_run_qs_check/
