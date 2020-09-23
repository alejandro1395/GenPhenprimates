#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/SBHs/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BBHs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/BlastP/
mkdir -p ${OUTDIR}

#create variable file
touch FinderBBH_input_JA.txt

#genes
#Loop for human genes
for filepath in $(ls ${INDIR}*/*tar.gz);
do echo "$filepath">>FinderBBH_input_JA.txt
echo $filepath
done
