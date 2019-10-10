#directory names of the exact & fuzz matching comparison sets
#e.g. fuzz0_no_a_post_check_flip
exact=$1
fuzz=$2

dir=$(dirname $0)

/bin/bash -x $dir/collect_all_stats.sh $exact
/bin/bash -x $dir/collect_all_stats.sh $fuzz

paste <(cut -d' ' -f 1,2,3 ${exact}.sorted.tsv) <(cut -d' ' -f 3 ${fuzz}.sorted.tsv) <(cut -d' ' -f 4 ${exact}.sorted.tsv) <(cut -d' ' -f 4 ${fuzz}.sorted.tsv) > both.sorted.pasted.tsv

# total_#multi_exonic_qisofrags        fuzz0_pct        fuzz20_pct        fuzz0_abs        fuzz20_abs
#(all_matching_contained_contains	all_matching  contained   contains    non_matching_overlaps   novel   repeats)
perl format_isoform_comparison_table.pl both.sorted.pasted.tsv  > both.sorted.pasted.formatted.tsv
