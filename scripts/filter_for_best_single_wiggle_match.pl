#!/bin/perl
#this takes in the partially processed and sorted output of the bedtools2 intersection
#of the source long reads junctions and a target set of junctions each with supporting # of reads
#input format (tab delimited):
#lr coords then target coords:
#chr    start   end num_reads   chr start   end target_name num_reads
#outputs the same format but only one line per unique LR junction

use strict;
use warnings;

my @cluster=();
while(my $line = <STDIN>)
{
    chomp($line); 
    my @f=split(/\t/,$line); 
    my ($c1,$s1,$e1,$nr1,$c2,$s2,$e2,$n2,$nr2)=@f;
    #is the current line overlapping the current cluster
    if(scalar(@cluster) > 0 && ($cluster[0]->[0] ne $c1 || $cluster[0]->[1] != $s1 ||  $cluster[0]->[2] != $e1))
    { 
        #now find the match with the >> # of read support (in the target)
        my $max = $cluster[0]; 
        for my $a (@cluster) 
        { 
            #track one with max # of reads
            $max=$a if($a->[8] > $max->[8]);
            #only trumped by one which exactly matches
            if($a->[1] == $a->[5] && $a->[2] == $a->[6])
            {
                $max=$a;
                last;
            }
        }
        print "".join("\t",@$max)."\n";
        @cluster=(); 
    }
    push(@cluster, \@f);
}

if(scalar(@cluster) > 0)
{
    my $max = $cluster[0]; 
    for my $a (@cluster) 
    { 
        #track one with max # of reads
        $max=$a if($a->[8] > $max->[8]);
        #only trumped by one which exactly matches
        if($a->[1] == $a->[5] && $a->[2] == $a->[6])
        {
            $max=$a;
            last;
        }
    }
    print "".join("\t",@$max)."\n";
}
