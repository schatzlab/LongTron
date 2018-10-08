#!/bin/perl
#like GenomeRanges.disjoin
use strict;

my $START=0;
my $END=1;

my $pc;
my $ps;
my $pe;
my @pa=();
my $payload="";
while(my $line=<STDIN>)
{
	chomp($line);
	my @f=split(/\t/,$line);
	my($c,$s,$e)=splice(@f,0,3);
	#overlap
	if($pc && $pc eq $c && $s <= $pe)
	{
		push(@pa,[$s,$START]);
		push(@pa,[$e,$END]);
		$pe = $e;
		$payload .= "\t".join("\t",@f);
		next;
	}
	if($pc)
	{
		disjoin($pc,$ps,$pe,\@pa,$payload);
	}
	$pc = $c;
	$ps = $s;
	$pe = $e;
	@pa = ([$s,$START],[$e,$END]);
	$payload = "\t".join("\t",@f);
}
if($pc)
{
	disjoin($pc,$ps,$pe,\@pa,$payload);
}

sub disjoin
{
	my ($c,$s,$e,$pa,$payload)=@_;
	my $pcoord;
	my $pt;
	my $i = 0;
	#sort by coordinate, this is both starts and ends
	#if coords are tied, favor STARTs over ENDs
	my @sa = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @$pa;
	my $len = scalar(@sa);
	for my $pos (@sa)
	{
		$i++;
		my ($coord, $t) = @$pos;

		#first in list
		my $range;
		if($pcoord)
		{
			#two starts in a row at the beginning are ignored (one is chosen)
			if($pt == $START && $t == $START && !($pcoord == $coord && $i == 2))
			{
				$range = "$pcoord\t".($coord-1);
			}
			elsif($pt == $START && $t == $END)
			{
				$range = "$pcoord\t$coord";
			}
			elsif($pt == $END && $t == $START)
			{
				$range = "".($pcoord+1)."\t".($coord-1);
			}
			#two ends in a row at the last are ignored (one is chosen)
			elsif($pt == $END && $t == $END && !($pcoord == $coord && $i == $len))
			{
				$range = "".($pcoord+1)."\t$coord";
			}
		}
		if($range)
		{
			print "$c\t$range".$payload."\n";
		}
		($pcoord, $pt) = ($coord, $t);
	}
}

sub disjoin_non_gr
{
	my ($c,$s,$e,$a,$payload)=@_;
	my $pcoord;
	my $i = 0;
	#sort by coordinate, this is both starts and ends
	my @sa = sort { $a <=> $b } @$a;
	my $len = scalar(@sa);
	for my $coord (@sa)
	{
		$i++;
		#avoid duplicating same coordinates
		#both current and the next one or the next and last one
		#but doesn't apply for a single, 1bp exon (s1,e1)
		next if($pcoord && $len > 2 && ($pcoord == $coord || ($i+1 == $len && $sa[$i] == $coord)));
		if($pcoord)
		{
			my $ccoord = ($i==$len?$coord:$coord-1);
			#off by one start for BED format
			print "$c\t".($pcoord-1)."\t$ccoord"."$payload\n";
		}
		$pcoord = $coord;
	}
}
