#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_out/

#We keep the species names for each one of the primates used for annotation from their path
Primates=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastN/
mkdir -p ${OUTDIR}Query/qu
mkdir -p ${OUTDIR}Query/out

for filepath in $(ls ${Primates}*/*.cds.*);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}Query/${species_name}/
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3
module load CDHIT/4.7

gunzip ${filepath};
cd-hit-est \
-i ${Primates}${species_name}/${species_name}.cds.fa \
-o ${OUTDIR}Query/${species_name}/${species_name}.nr100.cds.fa \
-c 1.0 \
-n 10 \
-M 16000 \
-d 0 \
-T 8;
gzip ${Primates}${species_name}/${species_name}.cds.fa;
gzip ${OUTDIR}Query/${species_name}/${species_name}.nr100.cds.fa" > ${OUTDIR}Query/qu/${species_name}.nr100.cds.sh
jobname=$(echo ${OUTDIR}Query/qu/${species_name}.nr100.cds.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}Query/out/${species_name}.nr100.cds.out \
-e ${OUTDIR}Query/out/${species_name}.nr100.cds.err -n ${species_name} -u 1 -t 1 -w 00:30:00
done
