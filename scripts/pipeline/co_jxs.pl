#!/usr/bin/perl
#find co-occurring jx/exons and print shared read (nro) counts as well as indiviudal (nr) counts
use strict;
use warnings;

my $jx2readsF = shift;

main();


sub map_reads_and_jxs($$$$)
{
	my ($line, $jx_id, $additional_args_ref) = @_;
	my ($jxs, $r2jx) = @$additional_args_ref;
	my ($c,$s,$e,$o,$nr,$reads) = split(/\t/,$line);
	my @reads = split(/,/,$reads);
	map { push(@{$r2jx->{$_}},$jx_id); } @reads;
	$jxs->{$jx_id} = [$c,$s,$e,$o,$nr];
}

sub map_jx_co($$$$$)
{
	my ($line, $jx_id, $additional_args_ref) = @_;
	my ($jxs, $r2jx, $print) = @$additional_args_ref;

	my ($c,$s,$e,$o,$nr,$reads) = split(/\t/,$line);
	my @reads = split(/,/,$reads);
	my %ojxs;
	#find every co-occurring jx from every read this jx occurs in
	for my $r (@reads)
	{
		my @ojxs = @{$r2jx->{$r}};
		for my $ojx (@ojxs)
		{
			next if($ojx == $jx_id);
			#track count of every shared read this other jx occurs in
			$ojxs{$ojx}++;
		}
	}
	if($print)
	{
		#now print mapping of this jx to all of its co-occuring jxs, sorted by the # of
		#reads it co-occurs in
		print "".scalar(keys %ojxs)."\t";
		print "".join("\t",($c,$s,$e,$o,$nr));
		#sort by # of reads each co-occurring jx appears in
		for my $ojx (sort { $ojxs{$a} <=> $ojxs{$b} } keys %ojxs)
		{
			my $nro = $ojxs{$ojx};
			my ($c2,$s2,$e2,$o2,$nr2,$nro2) = @{$jxs->{$ojx}};
			## # of shared co-occurrences, chrm, start, end, strand, # of supporting reads, 
			## # of co-occurrences for this co-occurring jx 
			print "\t".join("|",($nro,$c2,$s2,$e2,$o2,$nr2,$nro2));
		}
		print "\n";
		return;
	}
	#track number of co-occurring jxs for this jx
	push (@{$jxs->{$jx_id}}, scalar (keys %ojxs));
}

sub main
{
	my %r2jx;
	my %jxs;
	#first map:
	#1) read_id to jx_ids
	#2) jx_id to jx info
	process_input($jx2readsF, \&map_reads_and_jxs, \%jxs, \%r2jx);
	#now loop back through file
	# to actually get counts of all co-occurring jx for every jx
	process_input($jx2readsF, \&map_jx_co, \%jxs, \%r2jx);
	#this time print
	process_input($jx2readsF, \&map_jx_co, \%jxs, \%r2jx, 1);
}


sub process_input
{
	my ($filein, $process_fn) = (shift(@_), shift(@_));
	my $jx_id = 1;
	my $fin;
	open($fin,"<$filein") or die "$!\n";
	while(my $line=<$fin>)
	{
		chomp($line);
		$process_fn->($line, $jx_id, \@_);
		$jx_id++;
	}
	close($fin);
}
