#!/bin/perl
#for each long read, parse CIGAR to create a read-derived "isoform"
#where ends are assumed to be "fuzzy"
#primarily interested in order and location of junctions

use strict;
use warnings;

my $b = 0;
while(my $line = <STDIN>)
{
	chomp($line);
	#assume cut 1,2,3,4,6 has already been applied to BAM output
	my ($rname, $flag, $c, $s, $f) = split(/\t/,$line);
	my $o = (int($flag) & 0x10)?'-':'+';
	my $r = $s; 
	print "$rname\t$c\t$o\t$s\t";
	#print out list of junctions between start/end of "isoform"
	while($f=~/(\d+)([NMD=X])/cg) 
	{ 
		my $i=$1; 
		my $t=$2; 
		if($t eq "N") 
		{
			my $jx_s = $r;
			my $jx_e = $r+$i-1;
			print "$jx_s-$jx_e,";
		} 
		$r+=$i; 
	} 
	print "\t".($r-1)."\n";
	$b++;
}
