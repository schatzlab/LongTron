use strict;
use warnings;

my $refF = shift;
open(IN,"<$refF");
my %ref_introns;
my %ref_transcript_counts;
my ($pc,$ps,$pe,$po,$ptid);
while(my $line=<IN>)
{
    chomp($line);
    my @f=split(/\t/,$line);
    my ($c,$s,$e,$o,$tid)=($f[0],$f[3],$f[4],$f[6],$f[7]);
    #start of new transcript
    if(!$ptid || $ptid ne $tid)
    {
        $ptid=$tid;
        $pe=$e;
        next;
    }
    my $istart=$pe+1;
    my $iend=$s-1;
    my $key = join("|",($c,$istart,$iend,$o));
    push(@{$ref_introns{$key}},$tid);
    $ref_transcript_counts{$tid}++;
    $pe=$e;
}
close(IN);

my %tid_counts;
($pe,$ptid)=(undef,undef);
my $intron_count_per_query=0;
while(my $line = <STDIN>)
{
    chomp($line);
    my @f=split(/\t/,$line);
    my ($c,$s,$e,$o,$tid)=($f[0],$f[3],$f[4],$f[6],$f[7]);
    #start of new transcript
    if(!$ptid || $ptid ne $tid)
    {
        if($ptid && $intron_count_per_query > 0)
        {
            my $matches = 0;
            my $matching_tids = "";
            for my $rtid (keys %tid_counts)
            {
                #matching intron counts between the 2 chains have to match exactly
                if($tid_counts{$rtid} == $ref_transcript_counts{$rtid} && $tid_counts{$rtid} == $intron_count_per_query)
                {
                    $matches++;
                    $matching_tids .= ",$rtid";
                }
            }
            print "$ptid\t$matches\t$matching_tids\n";
            %tid_counts=();
        }
        $ptid=$tid;
        $pe=$e;
        $intron_count_per_query = 0;
        next;
    }
    $intron_count_per_query++;
    my $istart=$pe+1;
    my $iend=$s-1;
    my $key = join("|",($c,$istart,$iend,$o));
    my $matching_tids = $ref_introns{$key};
    map { $tid_counts{$_}++; } @$matching_tids;
}
     
if($ptid && $intron_count_per_query > 0)
{
    my $matches = 0;
    my $matching_tids = "";
    for my $rtid (keys %tid_counts)
    {
        #matching intron counts between the 2 chains have to match exactly
        if($tid_counts{$rtid} == $ref_transcript_counts{$rtid} && $tid_counts{$rtid} == $intron_count_per_query)
        {
            $matches++;
            $matching_tids .= ",$rtid";
        }
    }
    print "$ptid\t$matches\t$matching_tids\n";
    %tid_counts=();
}
