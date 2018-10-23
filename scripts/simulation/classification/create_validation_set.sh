#!/bin/bash
input=$1
cut -f 9- $input > ${input}.x
cut -f 8 ${input} | perl -ne 'BEGIN { %a=("problem-free"=>1,"non-recurrent"=>2,"recurrent"=>3,"novel"=>4); } chomp; $r=$_; @b=(0,0,0,0); $b[$a{$r}-1]=1; print "$r\t".join("\t",@b)."\n";' | cut -f 2- > ${input}.y
cut -f 8 ${input} | perl -ne 'chomp; $r=$_; if("problem_free" eq $r) { print "0\n"; next; } print "1\n";' > binary.y.validation
ln -fs ${input}.x validation.x
ln -fs ${input}.y validation.y
