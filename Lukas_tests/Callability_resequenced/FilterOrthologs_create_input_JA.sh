#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Callability_resequenced/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Lukas_tests/Callability_resequenced/

#create variable file
touch FilterOrthologs_create_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*.tsv);
do echo "$filepath">> FilterOrthologs_create_input_JA.txt
echo $filepath
done
