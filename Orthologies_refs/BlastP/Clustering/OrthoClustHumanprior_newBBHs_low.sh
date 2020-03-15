#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/

#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Primate_Traits/OUTPUT/TraitsPerSpecies.txt
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Orthology_clusters/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastP/Clustering/
mkdir -p ${OUTDIR}qu
mkdir -p ${OUTDIR}out
mkdir -p ${OUTDIR}OMAs

for filepath in $(ls ${INDIR}*/*.BBH.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${species_name}/
done
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}OrthoClustHumanprior_newBBHs_low.py ${INDIR} \
${SPECIES_IDs}summary_species.txt \
${TRAITS} \
${OUTDIR}" > ${OUTDIR}qu/all_species.pep.OrthoClustHumanprior_newBBHs_low.sh
jobname=$(echo ${OUTDIR}qu/all_species.pep.OrthoClustHumanprior_newBBHs_low.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/all_species.pep.OrthoClustHumanprior_newBBHs_low.out \
-e ${OUTDIR}out/all_species.pep.OrthoClustHumanprior_newBBHs_low.err -n Clustering -u 1 -t 1 -w 1-23:59:00 -r lowprio
