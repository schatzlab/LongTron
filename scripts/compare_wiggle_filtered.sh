#!/bin/bash
#parameters:
#1: min # of reads filter for LR
#2: file of matching junctions from pre-run matching w/ wiggle script (e.g. paftools.js), format: ^\s+#_reads_supporting_lr_jx\s+chr\tstart\tend\t#_reads_supporting_target_jx
#3: min # of reads filter for Target
#4: path to directory of already run exact, filtered matches (so we can get total wc numbers for filtered sets)

lr=`wc -l $4/lr.filtered.${1} | cut -d" " -f 1`
target=`wc -l $4/target.filtered.${3} | cut -d" " -f 1`

cat ${2} | perl -ne 'chomp; ($j,$lnr,$c,$s,$e,$tnr)=split(/\s+/,$_); if($lnr >= '${1}' && $tnr >= '${3}') { print "$c\t$s\t$e\n"; }' | wc -l | perl -ne 'chomp; $s=$_; $lr_per=$s/'${lr}'; $target_per=$s/'${target}'; print "'${1}'\t'${3}'\t$s\t'${lr}'\t'${target}'\t"; printf("%.3f\t%.3f\n",$lr_per,$target_per);'
