#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Prot_alignments_filtered/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Prot_alignments_qual/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_human_driven_refs/
mkdir -p ${OUTDIR}

#create variable file
touch trimAl_filter2_pos_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*/*.filter1.pep.aln);
do echo "$filepath">> trimAl_filter2_pos_input_JA.txt
echo $filepath
done
