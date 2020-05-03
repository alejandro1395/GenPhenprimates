#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Orthology_clusters/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/BlastP/
mkdir -p ${OUTDIR}

#create variable file
touch Clean_clusters_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*/*.gz);
do echo "$filepath">>Clean_clusters_input_JA.txt
echo $filepath
done
