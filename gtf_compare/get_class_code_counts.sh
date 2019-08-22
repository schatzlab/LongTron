find . -name "*.tmap" | perl -ne 'chomp; $s=$_; `wc -l $s | cut -d" " -f 1 > $s.wc`;'
find . -name "*.tmap" | perl -ne 'chomp; $s=$_; `cut -f 3 $s | sort | uniq -c > $s.counts`;'
