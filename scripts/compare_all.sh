#!/bin/bash

#for s in gtex short sra
for s in sra
do
	cd ${s}_exact_filtered
	#ln -fs ../${s}_junctions.tsv ./
	../run_exact_filtered.sh ../${s}_junctions.tsv > run2 2>&1
	cd ../${s}_exact_unfiltered
	../run_exact_unfiltered.sh ../${s}_exact_filtered/${s}_junctions.tsv.filtered.1 > run2 2>&1
	cd ../
done
s="refseq"
cd ${s}_exact_unfiltered
../run_exact_unfiltered.sh ../${s}_junctions.tsv > run2 2>&1
cd ../

