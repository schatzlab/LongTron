tmap_file=$1

#get all transcripts (including exons)
cut -f 3 $tmap_file | sort | uniq -c > ${tmap_file}.all_counts

cat ${tmap_file}.all_counts | perl -ne 'BEGIN { %m=("="=>"full_match","c"=>"match_in_cref","k"=>"match_in_kquery","p"=>"no_overlap","r"=>"no_overlap","u"=>"no_overlap"); } chomp; $f=$_; ($j,$cnt,$t)=split(/\s+/,$f); next if($t =~ /class_code/); $total+=$cnt; $category=($m{$t}?$m{$t}:"mismatching_overlap"); $h{$category}+=$cnt; END { @keys=sort { $a cmp $b } keys %h; print "$total\n"; print "".join(",",@keys)."\n"; print "".join(",",(map { $h{$_}; } @keys))."\n"; print "".join(",",(map { $v=$h{$_}; sprintf("%.1f\%",100*($v/$total)); } @keys))."\n";  }' > ${tmap_file}.all_counts.summarized

paste <(echo ${tmap_file}.all_counts | cut -d'/' -f 2) <(head -1 ${tmap_file}.all_counts.summarized) <(tail -n1 ${tmap_file}.all_counts.summarized)

#now get only intron chain (ic) counts
cut -f 3,6 $tmap_file | egrep -v -e '	1$' | cut -f 1 | sort | uniq -c > ${tmap_file}.ic_counts

cat ${tmap_file}.ic_counts | perl -ne 'BEGIN { %m=("="=>"full_match","c"=>"match_in_cref","k"=>"match_in_kquery","p"=>"no_overlap","r"=>"no_overlap","u"=>"no_overlap"); } chomp; $f=$_; ($j,$cnt,$t)=split(/\s+/,$f); next if($t =~ /class_code/); $total+=$cnt; $category=($m{$t}?$m{$t}:"mismatching_overlap"); $h{$category}+=$cnt; END { @keys=sort { $a cmp $b } keys %h; print "$total\n"; print "".join(",",@keys)."\n"; print "".join(",",(map { $h{$_}; } @keys))."\n"; print "".join(",",(map { $v=$h{$_}; sprintf("%.1f\%",100*($v/$total)); } @keys))."\n";  }' > ${tmap_file}.ic_counts.summarized

paste <(echo ${tmap_file}.ic_counts | cut -d'/' -f 2 | paste -d'.' - <(echo -n "ic")) <(head -1 ${tmap_file}.ic_counts.summarized) <(tail -n1 ${tmap_file}.ic_counts.summarized)
