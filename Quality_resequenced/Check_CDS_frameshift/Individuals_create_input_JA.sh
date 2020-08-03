#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/lkuderna/PGDP/GENE_VARIANTS_testings/combined/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_resequenced/Check_CDS_frameshift/

#create variable file
touch Individuals_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*tsv);
do echo "$filepath">> Individuals_input_JA.txt
echo $filepath
done

