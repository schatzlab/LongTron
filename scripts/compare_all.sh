#!/bin/bash
set -o pipefail -o nounset -o errexit 

for s in gtex short pacbio sra
do
	cd ${s}_exact_filtered
	ln -fs ../${s}_junctions.tsv ./
	../run_exact_filtered.sh ../${s}_junctions.tsv > run2 2>&1
	cd ../${s}_exact_unfiltered
	../run_exact_unfiltered.sh ../${s}_exact_filtered/target.filtered.1 > run2 2>&1
	cd ../
done
s="refseq"
cd ${s}_exact_unfiltered
../run_exact_unfiltered.sh ../${s}_junctions.tsv > run2 2>&1
cd ../

