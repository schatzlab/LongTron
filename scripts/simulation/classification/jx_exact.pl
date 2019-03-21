#!/usr/bin/env perl
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
    $small_count++;
} 
close(IN);
my $matching=0; 
while(my $z=<STDIN>)
{
    chomp($z); 
    my ($c,$s,$e)=split(/\t/,$z);
    my $k = join("\t",($c,$s,$e)); 
    #print "coord\t$z\t".join("\t",@{$h{$z}})."\n"; 
    if($h{$k})
    {
        $h{$k}->[3] += 1; 
        $matching++;
    }
}
#large_count will vary with the total number of matching recount jxs (it's not the total number of recount jxs)
print STDERR "".($matching/$small_count)."\t$matching\t$small_count\n";
my %h2;
#assign counts to read/transcript IDs
for my $a (keys %h) 
{ 
    my $c = $h{$a}->[3]; 
    if($c && $c > 0) 
    { 
        my @f1=split(/,/,$h{$a}->[2]); 
        for my $e (@f1) 
        { 
            my ($id)=split(/\./,$e);
            $h2{$id}+=1;
        } 
    } 
} 
#now print count for each read/transcript ID
#as well as ratio of splice matches to all splices in each TID
for my $tid (keys %h2)
{
    my $t_ratio = $h2{$tid}/$totals{$tid}; 
    print "$tid\t$t_ratio\t".$h2{$tid}."\t".$totals{$tid}."\n";
    print STDERR "$tid\t$t_ratio\t".$h2{$tid}."\t".$totals{$tid}."\n" if($t_ratio <= 0.5 && $totals{$tid} > 2);
}
