#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/trimAl/source/trimal
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Alignments_refs/Input_clusters/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Alignments_refs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_refs/Qualities/
mkdir -p ${OUTDIR}Prot_alignments_qual/summary_statistics/

for filepath in $(ls ${INDIR}*/*.clusters.pep.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}Prot_alignments_qual/summary_statistics/${species_name}
echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}ref_mean_ident_spc.py ${OUTDIR}Prot_alignments_qual/${species_name} \
${SPECIES_IDs}summary_species.txt \
${OUTDIR}Prot_alignments_qual/summary_statistics/${species_name}"> ${OUTDIR}Prot_alignments_qual/qu/${species_name}_mean_ident_spc.sh
jobname=$(echo ${OUTDIR}Prot_alignments_qual/qu/${species_name}_mean_ident_spc.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}Prot_alignments_qual/out/${species_name}_mean_ident_spc.out \
-e ${OUTDIR}Prot_alignments_qual/out/${species_name}_mean_ident_spc.err -n ${species_name} -u 1 -t 1 -w 03:00:00
done
