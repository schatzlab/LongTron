#e.g. /home/cwilks/gtf_compare/full_tmap_matching_fuzz20

#bdir=$(dirname $0)
compare_dir=$1
comparison=`echo "$compare_dir" | perl -ne 'chomp; $f=$_; @f=split(/\//,$f); $d=pop(@f);  print "$d\n";'`

#find $compare_dir -name "*.tmap" > ${comparison}.all_tmaps
find $compare_dir -name "comparison.tsv" > ${comparison}.all_comparisons
echo -n "" > ${comparison}.tsv
#for f in `cat ${comparison}.all_tmaps`; do /bin/bash ./make_isoform_comparison_table.sh $f >> ${comparison}.tsv ; done
find $compare_dir -name "comparison.tsv" -exec /bin/bash -c "cat <(echo -n {} | cut -d'/' -f 2 | tr '\n' ' ') <(fgrep -v 'Shell' {})" \; >> ${comparison}.tsv
fgrep "flipped" ${comparison}.tsv | sed -e 's/\.flipped//' > flipped
cut -d' ' -f 1 flipped > flipped.names
cat flipped <(fgrep -v -f flipped.names ${comparison}.tsv) | sort -k1,1 > ${comparison}.sorted.tsv

#paste <(cut -d' ' -f 1,2,3 ${comparison}.sorted.tsv) <(cut -d' ' -f 3 ${comparison}.sorted.tsv) <(cut -d' ' -f 4 ${comparison}.sorted.tsv) <(cut -d' ' -f 4 ${comparison}.sorted.tsv) > ${comparison}.sorted.pasted.tsv
