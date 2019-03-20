#!/usr/bin/env perl
use strict;
use warnings;

my $W = 20;

my $small_setF = shift;
my $windowed_mode = shift;

my $matching=0; 
my $large_count=0; 
my %h;
#e.g. trans_sim10.fl.junctions.clean.bed
open(IN,"<$small_setF");
my $small_count=0;
my %totals;
while(my $f=<IN>) 
{ 
    chomp($f); 
    my @f=(); 
    @f=split(/\t/,$f); 
    #don't adjust for BED, since both files will be BED
    #$f[1]++; 
    my $k=join("\t",splice(@f,0,3)); 
    #my @f1=splice(@f,0); 
    map { my ($id)=split(/\./,$_); $totals{$id}++; } split(/,/,$f[2]);
    $h{$k}=\@f;
    $small_count++;
} 
close(IN);
my $small_matches = 0;
my %matches;
while(my $z=<STDIN>)
{
    chomp($z); 
    $large_count++;
    if($h{$z} && !$windowed_mode) 
    { 
        #print "coord\t$z\t".join("\t",@{$h{$z}})."\n"; 
        $h{$z}->[3]++; 
        $matching++;
        next;
    }
    elsif($windowed_mode)
    { 
        #windowed mode
        my ($c,$s,$e,$o,$ct,$ids,$c2,$s2,$e2,$bps)=split(/\t/,$z);
        my $diff1 = abs($s2-$s);
        my $diff2 = abs($e-$e2);
        #both start/end need to be within 20 bases of overlapping start/end
        #and within the overlap (since we're running this with already widened windows)
        #if($bps > 0 && ($diff1 < $W && $diff1 >= 0) && ($diff2 < $W && $diff2 >= 0))
        if($bps > 0 && $diff1 < $W && $diff2 < $W)
        {
            $z = join("\t",($c,$s,$e)); 
            $matching++;
            #only count an overlapping match once
            if(!$matches{$z})
            {
                $small_matches++;
                $h{$z}->[3]++; 
            }
            $matches{$z}=1;
        }
    }
}
print STDERR "".($small_matches/$small_count)."\t".($matching/$large_count)."\t$small_matches\t$matching\t$small_count\t$large_count\n";
my %h2;
for my $a (keys %h) 
{ 
    my $c = $h{$a}->[3]; 
    if($c && $c > 0) 
    { 
        my @f1=split(/,/,$h{$a}->[2]); 
        for my $e (@f1) 
        { 
            my ($id)=split(/\./,$e);
            $h2{$id}+=$c;
        } 
    } 
} 
for my $tid (keys %h2)
{
    my $t_ratio = $h2{$tid}/$totals{$tid}; 
    print "$tid\t$t_ratio\t".$h2{$tid}."\t".$totals{$tid}."\n";
    print STDERR "$tid\t$t_ratio\t".$h2{$tid}."\t".$totals{$tid}."\n" if($t_ratio <= 0.5 && $totals{$tid} > 2);
}
