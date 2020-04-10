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
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_refs/Qualities/
mkdir -p ${OUTDIR}Prot_alignments_qual/
mkdir -p ${OUTDIR}Prot_alignments_qual/qu/
mkdir -p ${OUTDIR}Prot_alignments_qual/out/

echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}total_mean_ident_spc.py ${OUTDIR}Prot_alignments_qual \
${SPECIES_IDs}summary_species.txt \
${OUTDIR}Prot_alignments_qual/total_mean_ident_spc.tsv"> ${OUTDIR}Prot_alignments_qual/qu/total_mean_ident_spc.sh
jobname=$(echo ${OUTDIR}Prot_alignments_qual/qu/total_mean_ident_spc.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
#/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}Prot_alignments_qual/out/total_mean_ident_spc.out \
#-e ${OUTDIR}Prot_alignments_qual/out/total_mean_ident_spc.err -n total_mean_ident -u 1 -t 1 -w 06:00:00
