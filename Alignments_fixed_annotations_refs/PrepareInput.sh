#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/

#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Alignments_refs/Input_clusters/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_refs/
mkdir -p ${OUTDIR}Alignments_refs/Input_clusters/qu/
mkdir -p ${OUTDIR}Alignments_refs/Input_clusters/out/

for filepath in $(ls ${OUTDIR}*/*.clusters.pep.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}PrepareInput.py ${filepath} \
${SPECIES_IDs}summary_species.txt \
${SPECIES_IDs} \
${OUTDIR}${species_name}/${species_name}" > ${OUTDIR}qu/${species_name}.fastas.pep.sh
jobname=$(echo ${OUTDIR}qu/${species_name}.fastas.pep.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.fastas.pep.out \
-e ${OUTDIR}out/${species_name}.fastas.pep.err -n ${species_name} -u 1 -t 1 -w 23:59:00
done
