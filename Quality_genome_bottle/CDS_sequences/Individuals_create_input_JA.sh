#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/lkuderna/PGDP/project_NIST/GENE_VARIANTS/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_genome_bottle/CDS_sequences/

#create variable file
touch Individuals_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*/*private_CDS.tsv);
do echo "$filepath">> Individuals_input_JA.txt
echo $filepath
done

