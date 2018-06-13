#!/usr/bin/perl
use strict;
use warnings;

my $closest_outputF = shift;

#first read in closest/overlapping gene output
my %h; 
my %gene_strand; 
#open(IN,"<novel_jx_ranks.w_tcounts.sorted.closest");
open(IN,"<$closest_outputF");
while(my $line = <IN>) 
{
	chomp($line); 
	my @f=split(/\t/,$line);
	#ONT JX fields:
	#chrm, start, end, in_pacbio?, not_overlapping_exon?, num overlapping ONT reads
	my ($c,$s,$e,$pac,$ov,$nlr) = splice(@f,0,6);
	#adjust BED start coords
	$s++; 
	#closest/overlapping gene start
	$f[1]++; 
	#closet/overlapping gene distance (0 if overlapping)
	my $dist = pop(@f);
	#overwrite gene's chrm with gene ID since we already have the chrm
	$f[0] = $f[3]; 
	#remove old gene ID and unused score fields
	splice(@f,3,2); 
	#pack the of the gene fields into a string
	my $rest = join("|",@f); 
	#create JX coordinate hash key
	my $k = join(":",($c,$s,$e));
	#get gene strand	
	my $strand = $f[$#f]; 
	#multi gene
	map { $gene_strand{$k}->{$_}=1; } split(/,/,$strand);
	#new JX	
	if(!$h{$k})
	{
		$h{$k} = "".join("\t",($c,$s,$e,$pac,$ov,$nlr,$dist))."\t".$rest;
	}
	#just add this overlapping/closest gene to existing JX entry
	else
	{
		$h{$k} .= ",$rest";
	}
} 
close(IN);

#now match against all original ONT JXs to get additional fields
while(my $line=<STDIN>)
{
	chomp($line);
	my ($c,$s,$e,$o,$nr)=split(/\t/,$line);
	my $k=join(":",($c,$s,$e)); 
	my $gene_strands = join(",",sort { $a cmp $b } keys %{$gene_strand{$k}}); 
	my $f=$h{$k};
	#only interested in ones which are in the original list of ONT JXs
	if($f) 
	{ 
		my @f=split(/\t/,$f);
		#add in:
		#1) ONT jx strand
		#2) unique strand(s) of overlapping/closest genes
		#3) number of ONT reads *supporting* (not just overlapping) this JX
		$f[2].="\t$o\t$gene_strands\t$nr";
		print "".join("\t",@f)."\n"; 
	}
}
