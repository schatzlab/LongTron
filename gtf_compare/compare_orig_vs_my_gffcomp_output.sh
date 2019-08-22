#compare original (no changes) gffcompare vs. all mine changes including -a -T and --fuzzlength

#/path/to/orig_gffcompare_output/<dataset>.gtf.tmap
orig=$1
#/path/to/my_gffcompare_output/<dataset>.gtf.tmap
my=$2

export LC_ALL=C

cut -f 2,3,5,6 $orig | sort > orig.tmap.all.pairs
cut -f 1,3 orig.tmap.all.pairs | sort -u > orig.sorted.pairs
cut -f 2,3,5,6 $my | sort > mine.tmap.all.pairs
cut -f 1,3 mine.tmap.all.pairs | sort -u > mine.sorted.pairs

comm -1 -2 orig.sorted.pairs mine.sorted.pairs > sorted.pairs.shared
comm -2 -3 orig.sorted.pairs mine.sorted.pairs > sorted.pairs.only_in.orig
comm -1 -3 orig.sorted.pairs mine.sorted.pairs > sorted.pairs.only_in.mine

#make set of query ids unique by adding tabs
cut -f 2 sorted.pairs.shared | sed -e 's/$/\t/' | sed -e 's/^/\t/' > sorted.pairs.shared.q
cut -f 2 sorted.pairs.only_in.orig | sed -e 's/$/\t/' | sed -e 's/^/\t/' > sorted.pairs.only_in.orig.q
cut -f 2 sorted.pairs.only_in.mine | sed -e 's/$/\t/' | sed -e 's/^/\t/' > sorted.pairs.only_in.mine.q

#this should be most, the few which aren't in "mine" but in "orig" will be ones that were dup filtered out by the symmetric filtering (-a) I've added
#which was originally only done on the reference side. A wayt to check this is to re-run gffcompare (orig) with the two datasets reversed, then 
#see if any of the queries which are here only orig are in references list of tmap file (none should be)
fgrep -f sorted.pairs.only_in.orig.q mine.tmap.all.pairs > sorted.pairs.only_in.orig.q.in_mine
#all of these should be in original
fgrep -f sorted.pairs.only_in.mine.q orig.tmap.all.pairs > sorted.pairs.only_in.mine.q.in_orig

#do the main category change comparison on the shared pairs (should be vast majority)
fgrep -f sorted.pairs.shared.q orig.tmap.all.pairs | sort > orig.tmap.all.pairs.shared
fgrep -f sorted.pairs.shared.q mine.tmap.all.pairs | sort > mine.tmap.all.pairs.shared

comm -1 -2 orig.tmap.all.pairs.shared mine.tmap.all.pairs.shared > tmap.all.pairs.shared.same
comm -2 -3 orig.tmap.all.pairs.shared mine.tmap.all.pairs.shared > tmap.all.pairs.shared.orig
comm -1 -3 orig.tmap.all.pairs.shared mine.tmap.all.pairs.shared > tmap.all.pairs.shared.mine

paste tmap.all.pairs.shared.orig tmap.all.pairs.shared.mine > tmap.all.pairs.shared.diff.pasted
#the vast majority of the category changes should be an improvement orig=>mine meaning a better category (e.g. "j"->"=")
cut -f 2,6 tmap.all.pairs.shared.diff.pasted | sort | uniq -c > tmap.all.pairs.shared.diff.pasted.cat_counts
