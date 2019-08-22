count_file=$1

cat $count_file | perl -ne 'BEGIN { %m=("="=>"full_match","c"=>"match_in_cref","k"=>"match_in_kquery","p"=>"no_overlap","r"=>"no_overlap","u"=>"no_overlap"); } chomp; $f=$_; ($j,$cnt,$t)=split(/\s+/,$f); next if($t =~ /class_code/); $total+=$cnt; $category=($m{$t}?$m{$t}:"mismatching_overlap"); $h{$category}+=$cnt; END { @keys=sort { $a cmp $b } keys %h; print "$total\n"; print "".join(",",@keys)."\n"; print "".join(",",(map { $h{$_}; } @keys))."\n"; print "".join(",",(map { $v=$h{$_}; sprintf("%.1f\%",100*($v/$total)); } @keys))."\n";  }' > ${count_file}.summarized


paste <(echo $count_file | cut -d'/' -f 1) <(head -1 ${count_file}.summarized) <(tail -n1 ${count_file}.summarized)
