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
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastP/
mkdir -p ${OUTDIR}Query/qu
mkdir -p ${OUTDIR}Query/out

for filepath in $(ls ${Primates}*/*.pep.*);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}Query/${species_name}/
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}Filternrpep.py ${filepath} ${OUTDIR}Query/${species_name}/${species_name}.nr.pep.fa.gz" > ${OUTDIR}qu/${species_name}.nr.pep.sh
jobname=$(echo ${OUTDIR}qu/${species_name}.nr.pep.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.nr.pep.out \
-e ${OUTDIR}out/${species_name}.nr.pep.err -n ${species_name} -u 1 -t 1 -w 00:05:00
done
