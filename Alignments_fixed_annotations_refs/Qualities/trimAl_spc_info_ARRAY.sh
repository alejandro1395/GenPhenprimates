#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/trimAl/source/trimal
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Alignments_refs/Prot_alignments/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Alignments_refs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_refs/


# grab the files, and export it so the ‘child’ sbatch jobs can access it
export FILES=($(ls -1 ${OUTDIR}Prot_alignments_qual/qu/*.trimAl_spc_info.sh))
#get size of array
NUMFILES=${#FILES[@]}
#change first index
ZBNUMFILES=$(($NUMFILES – 1))
# now submit to SLURM
if [ $ZBNUMFILES -ge 0 ]; then
sbatch –array=0-$ZBNUMFILES trimAl_spc_info_ARRAY.sbatch
fi 
