#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/trimAl/source/trimal
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/Prot_alignments_filtered/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/Prot_alignments_qual/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_fixed_annotations_refs/Qualities/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Traits_counts_refs.tsv
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/summary_statistics/
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/summary_statistics/out/
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/summary_statistics/qu/
SUMMARY=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/summary_statistics/

echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}total_mean_ident_spc.py ${OUTDIR} \
${SPECIES_IDs}summary_species.txt \
${SUMMARY}"> ${SUMMARY}qu/total_mean_ident_spc.sh
jobname=$(echo ${SUMMARY}qu/total_mean_ident_spc.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${SUMMARY}out/total_mean_ident_spc.out \
-e ${SUMMARY}out/total_mean_ident_spc.err -n total_mean_ident -u 1 -t 1 -w 03:00:00
