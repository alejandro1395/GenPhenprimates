#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Prot_alignments_filtered/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Input_fastas/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Genes_predictors/PGLS/

#create variable file
touch SortCDS_repair_JA.txt

#Loop for human genes
genes="ART5 CNN3 CREBZF FBXO27 HSPB11 IQCF6 MDK MYL10 OR10AC1 OR5H8 "
echo $genes | tr " " "\n" | while read gene;
do for filepath in $(ls ${INDIR}${gene}/*.tar.gz);
do echo "$filepath">> SortCDS_repair_JA.txt
echo $filepath
done
done
