#!/usr/bin/perl
#for each long read, parse CIGAR to create a read-derived "isoform"
#where ends are assumed to be "fuzzy"
#primarily interested in order and location of junctions

use strict;
use warnings;

my $DEFAULT_CUFF_INFO='FPKM 100.0; frac 1.0; conf_lo 0.01; conf_hi 0.98; cov 101.0; full_read_support "yes";';

my $rid = 0;
while(my $line = <STDIN>)
{
	chomp($line);
	$rid++;
	#assume cut 1,2,3,4,6 has already been applied to BAM output
	my ($rname, $flag, $c, $s, $f) = split(/\t/,$line);
	my $o = (int($flag) & 0x10)?'-':'+';
	my $ps = $s;
	my $e = $s; 
	print "$c\tCufflinks\ttranscript\t$s";
	my $tlength = 0;
	my $exons = "";
	#print out list of junctions between start/end of "isoform"
	while($f=~/(\d+)([NMD=X])/cg) 
	{ 
		my $i=$1; 
		$tlength += $i;
		my $t=$2; 
		if($t eq "N") 
		{
			#exon coordinates
			$exons .= "$c\tCufflinks\texon\t$ps\t".($e-1)."\t1000\t$o\t.\tgene_id \"$rid\"; transcript_id \"$rid.1\"; $DEFAULT_CUFF_INFO\n";
			#junction coordinates
			#my $jx_s = $e;
			#my $jx_e = $e+$i-1;
			#print "$jx_s-$jx_e,";
			$e+=$i; 
			$ps = $e;
		}
		else
		{
			$e+=$i;
		}
	} 
	$exons .= "$c\tCufflinks\texon\t$ps\t".($e-1)."\t1000\t$o\t.\tgene_id \"$rid\"; transcript_id \"$rid.1\"; $DEFAULT_CUFF_INFO\n";

	print "\t".($e-1)."\t1000\t$o\t.\tgene_id \"$rid\"; transcript_id \"$rid.1\"; $DEFAULT_CUFF_INFO\n";
	print "$exons";
	#print "\t".($e-1)."\n";
}
