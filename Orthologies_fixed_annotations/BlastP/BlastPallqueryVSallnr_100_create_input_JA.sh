#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
QUERY=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_in/
TARGET=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_nrDB/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_out/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_fixed_annotations/BlastP/
mkdir -p ${OUTDIR}BLAST_out

#create variable file
touch BlastPallqueryVSallnr_100_input_JA.txt

#Loop for human genes
for filepath in $(ls ${QUERY}*/*.gz);
do echo "$filepath">>BlastPallqueryVSallnr_100_input_JA.txt
echo $filepath
mkdir -p ${OUTDIR}BLAST_nrDB/${gene_name}
done
