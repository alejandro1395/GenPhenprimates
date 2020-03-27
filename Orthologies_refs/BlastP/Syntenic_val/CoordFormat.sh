#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/

#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
CHROM_path=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Syntenic_alignments
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Orthology_synteny_val/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastP/Syntenic_val/
mkdir -p ${OUTDIR}qu
mkdir -p ${OUTDIR}out

for filepath in $(ls ${OUTDIR}*/*Validation_clusters.pep.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${species_name}/
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}CoordFormat.py ${filepath} \
${SPECIES_IDs}summary_species.txt \
${SPECIES_IDs} \
${CHROM_path} \
${OUTDIR}${species_name}/${species_name}.ref." > ${OUTDIR}qu/${species_name}.Clusters_loc.pep.sh
jobname=$(echo ${OUTDIR}qu/${species_name}.Clusters_loc.pep.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.Clusters_loc.pep.out \
-e ${OUTDIR}out/${species_name}.Clusters_loc.pep.err -n ${species_name} -u 1 -t 1 -w 06:00:00
done
