#!/usr/bin/perl
#find co-occurring jx/exons and print shared read (nro) counts as well as indiviudal (nr) counts
use strict;
use warnings;

my $jx2readsF = shift;

main();

sub main
{
	open(IN,"<$jx2readsF");
	my %r2jx;
	my %jxs;
	my $jx_id = 1;
	#first map:
	#1) read_id to jx_ids
	#2) jx_id to jx info
	while(my $line=<IN>)
	{
		chomp($line);
		my ($c,$s,$e,$o,$nr,$reads) = split(/\t/,$line);
		my @reads = split(/,/,$reads);
		map { push(@{$r2jx{$_}},$jx_id); } @reads;
		$jxs{$jx_id}=[$c,$s,$e,$o,$nr];
		$jx_id++;
	}
	close(IN);
	# to actually get counts of all co-occurring jx for every jx
	open(IN,"<$jx2readsF");
	#now loop back through file
	$jx_id = 1;
	while(my $line=<IN>)
	{
		chomp($line);
		my ($c,$s,$e,$o,$nr,$reads) = split(/\t/,$line);
		my @reads = split(/,/,$reads);
		my %ojxs;
		#find every co-occurring jx from every read this jx occurs in
		for my $r (@reads)
		{
			my @ojxs = @{$r2jx{$r}};
			for my $ojx (@ojxs)
			{
				next if($ojx == $jx_id);
				#track count of every shared read this other jx occurs in
				$ojxs{$ojx}++;
			}
		}
		#now print mapping of this jx to all of its co-occuring jxs, sorted by the # of
		#reads it co-occurs in
		print "".scalar(keys %ojxs)."\t";
		print "".join("\t",($c,$s,$e,$o,$nr));
		#sort by # of reads each co-occurring jx appears in
		for my $ojx (sort { $ojxs{$a} <=> $ojxs{$b} } keys %ojxs)
		{
			my $nro = $ojxs{$ojx};
			my ($c2,$s2,$e2,$o2,$nr2) = @{$jxs{$ojx}};
			print "\t".join("|",($nro,$c2,$s2,$e2,$o2,$nr2));
		}
		print "\n";
		$jx_id++;
	}
	close(IN);
}
