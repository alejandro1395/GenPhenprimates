#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/trimAl/source/trimal
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Alignments_refs/Prot_alignments/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_refs/
mkdir -p ${OUTDIR}qu/
mkdir -p ${OUTDIR}out/

for filepath in $(ls ${OUTDIR}*/*.pep.aln);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
gene_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 2)
mkdir -p ${OUTDIR}Prot_alignments/${species_name}/
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

${BIN} -in ${filepath} \
-sident > ${OUTDIR}${species_name}/${species_name}.${gene_name}.SpeciesIdentities" > ${OUTDIR}qu/${species_name}.${gene_name}.trimAl_spc_info.sh
jobname=$(echo ${OUTDIR}qu/${species_name}.${gene_name}.trimAl_spc_info.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
#/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.${gene_name}.trimAl_spc_info.out \
#-e ${OUTDIR}out/${species_name}.${gene_name}.trimAl_spc_info.err -n ${species_name}.${gene_name} -u 1 -t 1 -w 00:05:00
done
