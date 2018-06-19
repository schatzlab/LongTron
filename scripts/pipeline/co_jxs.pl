#!/usr/bin/perl
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
	open(IN,"<$jx2readsF");
	$jx_id = 1;
	while(my $line=<IN>)
	{
		chomp($line);
		my ($c,$s,$e,$o,$nr,$reads) = split(/\t/,$line);
		my @reads = split(/,/,$reads);
		my %ojxs;
		for my $r (@reads)
		{
			my @ojxs = @{$r2jx{$r}};
			for my $ojx (@ojxs)
			{
				next if($ojx == $jx_id);
				$ojxs{$ojx}++;
			}
		}
		#now print mapping of this jx to all of its co-occuring jxs, sorted by the number of
		#reads it co-occurs in
		print "".scalar(keys %ojxs)."\t";
		print "".join("\t",($c,$s,$e,$o,$nr));
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



		

