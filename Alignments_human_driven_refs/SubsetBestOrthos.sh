#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/

#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Orthology_clusters/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/
mkdir -p ${OUTDIR}Alignments_refs/
mkdir -p ${OUTDIR}Alignments_refs/Input_clusters/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_refs/
mkdir -p ${OUTDIR}Alignments_refs/Input_clusters/qu/
mkdir -p ${OUTDIR}Alignments_refs/Input_clusters/out/

for filepath in $(ls ${INDIR}*/*clean.pep.fast_fast.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}Alignments_refs/Input_clusters/${species_name}/
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}SubsetBestOrthos.py ${filepath} \
${SPECIES_IDs}summary_species.txt \
${OUTDIR}Alignments_refs/Input_clusters/${species_name}/${species_name}.clusters.pep.gz" > ${OUTDIR}Alignments_refs/Input_clusters/qu/${species_name}.clusters.pep.gz.sh
jobname=$(echo ${OUTDIR}Alignments_refs/Input_clusters/qu/${species_name}.clusters.pep.gz.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}Alignments_refs/Input_clusters/out/${species_name}.clusters.pep.gz.out \
-e ${OUTDIR}Alignments_refs/Input_clusters/out/${species_name}.clusters.pep.gz.err -n ${species_name} -u 1 -t 1 -w 06:00:00
done
