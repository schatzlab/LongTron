#e.g. /home/cwilks/gtf_compare/full_tmap_matching_fuzz20

#bdir=$(dirname $0)
compare_dir=$1
comparison=`echo "$compare_dir" | perl -e 'chomp; $f=$_; @f=split(/\//,$f); $d=pop(@f);  print "$d\n";'`

find $compare_dir -name "*.tmap" > ${comparison}.all_tmaps
echo -n "" > ${comparison}.tsv
for f in `cat ${comparison}.all_tmaps`; do /bin/bash ./make_isoform_comparison_table.sh $f >> ${comparison}.tsv ; done
sort -k1,1 ${comparison}.tsv > ${comparison}.sorted.tsv
