#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS/CDSsequences/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Input_fastas/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Genes_predictors/PGLS/

#create variable file
touch AlignCodons_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*/*.sorted.cds.fa.gz);
do echo "$filepath">> AlignCodons_input_JA.txt
echo $filepath
done
