#!/bin/bash
set -o pipefail -o nounset -o errexit 

#assumes samtools is in the path
MINIMAP=minimap2
SURVIVOR=./SURVIVOR/Debug/SURVIVOR
SURVIVOR_FL=./SURVIVOR_fl/Debug/SURVIVOR
GFFREAD=./gffread/gffread
HG38=/scratch/groups/blangme2/indexes/minimap2/hg38_plus.fa

#create FASTA of GencodeV28 annotation
${GFFREAD} -w gencode.v28.basic.annotation.fa -g GRCh38_full_analysis_set_plus_decoy_hla.fa gencode.v28.basic.annotation.gtf
#map ONT DirectRNA FASTQ against GencodeV28 annotation FASTA to help create transcript error profile
${MINIMAP2} --MD -L -ax map-ont gencode.v28.basic.annotation.fa NA12878-DirectRNA.pass.dedup.fastq.gz | samtools view -b - > NA12878-DirectRNA.vs.gencode_v28_basic.bam
#create error profile with minimum length 100 read
samtools view NA12878-DirectRNA.vs.gencode_v28_basic.sorted.bam | ${SURVIVOR} scanreads 100 survivor_errors100.all.txt
#generate simulated transcripts with DirectRNA error profile 10x coverage
${SURVIVOR} simreads gencode.v28.basic.annotation.fa survivor_errors100.all.txt 10 hg38_gencode_based_transcripts_sim10.fa
#align simulated transcripts against HG38
${MINIMAP2} -ax splice -uf -k14 -t 16 ${HG38} hg38_gencode_based_transcripts_sim10.fa | samtools view -b - | samtools sort -O bam -T  hg38_gencode_based_transcripts_sim10.1 - | tee hg38_gencode_based_transcripts_sim10.fa.bam | samtools index
mv -- -.bai hg38_gencode_based_transcripts_sim10.fa.bam.bai
/bin/bash -x extract.sh hg38_gencode_based_transcripts_sim10.fa.bam gencode.v28.basic.annotation.gtf

#alternately generate ~full length simulated transcripts
${SURVIVOR_FL} simreads gencode.v28.basic.annotation.fa survivor_errors100.all.txt 10 hg38_gencode_based_transcripts_sim10.fl.fa
${MINIMAP2} -ax splice -uf -k14 -t 16 ${HG38} hg38_gencode_based_transcripts_sim10.fa.fl | samtools view -b - | samtools sort -O bam -T  hg38_gencode_based_transcripts_sim10.fl.1 - | tee hg38_gencode_based_transcripts_sim10.fl.fa.bam | samtools index
mv -- -.bai hg38_gencode_based_transcripts_sim10.fl.fa.bam.bai

#get exons, raw-read based isoforms, do cuff compare on them, and create BigBed for visualization of raw isoforms
/bin/bash -x extract.sh hg38_gencode_based_transcripts_sim10.fl.fa.bam gencode.v28.basic.annotation.gtf
