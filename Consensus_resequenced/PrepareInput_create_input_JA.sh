#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/RESEQUENCED/Input_fastas/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Consensus_resequenced/
mkdir -p ${OUTDIR}

#create variable file
touch PrepareInput_create_input_JA.txt

#Loop for human genes
for filepath in $(ls ${INDIR}*/*.fasta);
do echo "$filepath">> PrepareInput_input_JA.txt
echo $filepath
done
