#!/usr/bin/env perl
#takes a full list of snaptron or annotated junctions and matches them using exact coordinate matches
#to the set passed in as a file on the command line
use strict;
use warnings;

my $small_setF = shift;

my %h;
#e.g. trans_sim10.fl.junctions.clean.bed
open(IN,"<$small_setF");
my $small_count=0;
my %totals;
while(my $f=<IN>) 
{ 
    chomp($f); 
    my @f=split(/\t/,$f); 
    #don't adjust for BED, since both files will be BED
    #$f[1]++; 
    my $k=join("\t",splice(@f,0,3)); 
    #my @f1=splice(@f,0); 
    map { my ($id)=split(/\./,$_); $totals{$id}++; } split(/,/,$f[2]);
    $h{$k}=\@f;
    $h{$k}->[3]=0;
    $h{$k}->[4]=0;
    $small_count++;
} 
close(IN);
my $matching=0; 
#read snaptron/annotation
while(my $z=<STDIN>)
{
    chomp($z); 
    my ($c,$s,$e,$o,$a)=split(/\t/,$z);
    my $k = join("\t",($c,$s,$e)); 
    #print "coord\t$z\t".join("\t",@{$h{$z}})."\n"; 
    if($h{$k})
    {
        $h{$k}->[3] += 1;
        $h{$k}->[4] = 1 if($a != 0);
    }
}
my %h2;
my $num_annotated = 0;
#assign counts to read/transcript IDs
for my $k (keys %h) 
{ 
    my $c = $h{$k}->[3]; 
    my $a = $h{$k}->[4];
    if($c && $c > 0) 
    { 
        $matching++;
        $num_annotated++ if($a != 0);
        my @f1=split(/,/,$h{$k}->[2]); 
        for my $e (@f1) 
        { 
            my ($id)=split(/\./,$e);
            $h2{$id}+=1;
        } 
    } 
} 
my $novel = $small_count - $matching;
my $snaptron_unannotated = $matching - $num_annotated;
#printf("%.2f\t%d\t%d\t%.2f\t%d\t%d\t%.2f\t%d\t%d\n",($novel/$small_count),$novel,$small_count,($num_annotated/$small_count),$num_annotated,$small_count,($snaptron_unannotated/$small_count),$snaptron_unannotated,$small_count);
printf("%d (100\%)\t%d (%.0f\%)\t%d (%.0f\%)\t%d (%.0f\%)\n",$small_count,$num_annotated,100*($num_annotated/$small_count),$matching,100*($matching/$small_count),$novel,100*($novel/$small_count));

#now print count for each read/transcript ID
#as well as ratio of splice matches to all splices in each TID
for my $tid (keys %h2)
{
    my $t_ratio = $h2{$tid}/$totals{$tid}; 
    print "$tid\t$t_ratio\t".$h2{$tid}."\t".$totals{$tid}."\n";
    print STDERR "$tid\t$t_ratio\t".$h2{$tid}."\t".$totals{$tid}."\n" if($t_ratio <= 0.5 && $totals{$tid} > 2);
}
