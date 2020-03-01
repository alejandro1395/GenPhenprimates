#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/

#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastP/BBHs/
mkdir -p ${OUTDIR}qu
mkdir -p ${OUTDIR}out

for filepath in $(ls ${OUTDIR}*/*.SBH.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${species_name}/
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}FinderBBH.py ${OUTDIR} ${species_name} \
${SPECIES_IDs}summary_species.txt \
${OUTDIR}${species_name}/${species_name}.pep.BBH.gz" > ${OUTDIR}qu/${species_name}.pep.BBH.sh
jobname=$(echo ${OUTDIR}qu/${species_name}.pep.BBH.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.pep.BBH.out \
-e ${OUTDIR}out/${species_name}.pep.BBH.err -n ${species_name} -u 1 -t 1 -w 06:00:00
done
