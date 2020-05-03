#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Input_fastas/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Prot_alignments/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_human_driven_refs/
mkdir -p ${OUTDIR}

#create variable file
touch MSA_create_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*/*.gz);
do echo "$filepath">> MSA_create_input_JA.txt
echo $filepath
done
