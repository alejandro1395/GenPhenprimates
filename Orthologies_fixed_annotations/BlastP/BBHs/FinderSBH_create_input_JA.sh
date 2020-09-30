#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_out/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/SBHs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_fixed_annotations/BlastP/BBHs/
mkdir -p ${OUTDIR}

#create variable file
touch FinderSBH_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*/*tar.gz);
do echo "$filepath">>FinderSBH_input_JA.txt
echo $filepath
done
