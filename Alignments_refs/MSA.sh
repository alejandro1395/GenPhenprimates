#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/MUSCLE/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Alignments_refs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_refs/
mkdir -p ${OUTDIR}Prot_alignments/
mkdir -p ${OUTDIR}Prot_alignments/qu/
mkdir -p ${OUTDIR}Prot_alignments/out/

for filepath in $(ls ${OUTDIR}Input_clusters/*/*.pep.fa.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
gene_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 2)
mkdir -p ${OUTDIR}Prot_alignments/${species_name}/
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

gunzip -c ${filepath} | ${BIN}muscle3.8.31_i86linux64 -in - \
-clwout ${OUTDIR}Prot_alignments/${species_name}/${species_name}.${gene_name}.pep.aln" > ${OUTDIR}Prot_alignments/qu/${species_name}.${gene_name}.pep.sh
jobname=$(echo ${OUTDIR}Prot_alignments/qu/${species_name}.${gene_name}.pep.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}Prot_alignments/out/${species_name}.${gene_name}.pep.out \
-e ${OUTDIR}Prot_alignments/out/${species_name}.${gene_name}.pep.err -n ${species_name}.${gene_name} -u 1 -t 1 -w 00:05:00
done
