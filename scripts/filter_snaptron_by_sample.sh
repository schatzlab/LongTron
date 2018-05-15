#!/bin/bash
#takes the rail sample id (rail_id) of the sample you want to filter for in the Snaptron data
#also takes the minimum read coverage for the sample for that junction

cat /dev/stdin | perl -ne 'chomp; $f=$_; @f=split(/\t/,$f); if($f[11]=~/,'${1}':(\d+)/) { $n=$1; next if($n < '${2}'); print "$f\n";}'
