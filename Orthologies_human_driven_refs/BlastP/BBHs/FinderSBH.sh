#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/

#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
INDIR_CLUST=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/nrDB/
INDIR_BLAST=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_out/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastP/BBHs/
mkdir -p ${OUTDIR}qu
mkdir -p ${OUTDIR}out

for filepath in $(ls ${INDIR_BLAST}*/*.pep.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${species_name}/
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}FinderSBH.py ${filepath} ${INDIR_CLUST}allspeciesDBnr_100.pep.fa.clstr \
${SPECIES_IDs}summary_species.txt \
${OUTDIR}${species_name}/${species_name}.pep.SBH.gz" > ${OUTDIR}qu/${species_name}.pep.SBH.sh
jobname=$(echo ${OUTDIR}qu/${species_name}.pep.SBH.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.pep.SBH.out \
-e ${OUTDIR}out/${species_name}.pep.SBH.err -n ${species_name} -u 1 -t 1 -w 06:00:00
done
