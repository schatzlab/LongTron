
#e.g. 0 for all_matching_contained_contains %
my $PCT_COL=0;

my %col_map = ("na12878_G029"=>0, "na12878_all"=>1, "skbr3_G029"=>0, "skbr3_all"=>1, "na12878_short"=>2, "skbr3_short"=>2, "na12878_ox"=>3,"na12878_oxflair"=>3,"na12878_oxraw"=>3,"na12878_pb"=>3,"na12878_pbflair"=>3,"na12878_pbraw"=>3,"skbr3_pb"=>3);

my %name_map = ("na12878_ox"=>["Oxford NA12878 (OX-RAW)", "OX-RAW"], "na12878_oxraw"=>["Oxford NA12878 (OX-RAW)", "OX-RAW"], "na12878_pb"=>["PacBio NA12878 (PB-RAW)", "PB-RAW"], "na12878_pbraw"=>["PacBio NA12878 (PB-RAW)", "PB-RAW"], "skbr3_pb"=>["PacBio SKBR3 (PB-SKBR3)", "PB-SKBR3"], "na12878_short"=>["Illumina NA12878", "Illumina"], "skbr3_short"=>["Illumina SKBR3", "Illumina"], "na12878_pbflair"=>["PacBio NA12878 FLAIR", "PB-FLAIR"], "na12878_oxflair"=>["Oxford NA12878 FLAIR", "OX-FLAIR"]);
 

#e.g. fuzz20_no_a_post_check_flip.both.sorted.pasted.tsv
my $filein=shift;
open(IN,"<$filein");

my %totals;

print "Sample\tTotal Intron Chains\tGencode V29\tUnion of Annotations\tShort Assembly\t(Other) Long Reads\n";

while(my $line=<IN>)
{
    chomp($line);
    my ($name, $total_isos, $exact_pct, $fuzz_pct, $exact_abs, $fuzz_abs) = split(/\s+/,$line);
    my @name = split(/_/,$name);
    if($name =~ /_((all)|(G029))_/ || scalar (@name) == 3)
    {
        #TODO special handling of annotation only
        next;
    }
    my ($tech, $query_set, $vs, $ref_set) = @name;
    my $qname = $tech."_".$query_set;
    my $rname = $tech."_".$ref_set;
    
    my ($qfull, $qshort) = @{$name_map{$qname}};
    my ($rfull, $rshort) = @{$name_map{$rname}};
    my $rcol = $col_map{$rname};
   
    #my ($all_matching_contained_contains, $all_matching, $contained, $contains, $non_matching_overlaps, $novel, $repeats) = split(/,/,$exact_pct);
    my @exact_pct = split(/,/,$exact_pct);
    my $exact_val = $exact_pct[$PCT_COL];
    my @fuzz_pct = split(/,/,$fuzz_pct);
    my $fuzz_val = $fuzz_pct[$PCT_COL];
    #TODO: keep total_isos for each set
  
    my $value = "$rshort: $exact_val ($fuzz_val)"; 
    $value = "$exact_val ($fuzz_val)" if($ref_set =~ /(all)|(G029)/);
    push(@{$totals{$qfull}->[$rcol]}, $value);
}
close(IN);

for my $k (sort { $b cmp $a } keys %totals)
{
    my $out="$k\t";
    my $v = $totals{$k};
    for my $vs (@$v)
    {
        $vs = ["NA"] if(scalar(@$vs) == 0);
        $out .= "".join(", ",@$vs)."\t";
    }
    $out =~ s/\t$//;
    $out .= "\tNA" if($k eq $name_map{'skbr3_pb'}->[0]);
    print "$out\n";
}
