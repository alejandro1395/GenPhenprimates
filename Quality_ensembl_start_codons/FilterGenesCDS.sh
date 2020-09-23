#!/bin/bash
#SBATCH --array=1-2
#SBATCH --job-name=QualCDS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=23:59:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BEDTOOLS

#Define PATH argument
SPECIES_GFF=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/GFFs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_ensembl_start_codons/
HUMAN_MODELS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/Human_models/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/GFFs/
REF_FASTA=/scratch/devel/lkuderna//PGDP/REFERENCES_FROZEN/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}GFFs_species.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
zcat ${SPECIES_GFF}${ARGUMENT1} | awk -v OFS='\t' -v FS='\t' '$3 == "CDS" {print $0}' > \
${OUTDIR}${species_name}.CDS.gff;
gzip ${OUTDIR}/${species_name}.CDS.gff;
#IMPORT ZIPPED FILE INTO BEDTOOLS GETFASTA
zcat ${OUTDIR}/${species_name}.CDS.gff.gz | bedtools getfasta -fi ${REF_FASTA}${species_name}/${species_name}.fasta/*.fasta -bed - -fo stdout | \
paste - - > ${OUTDIR}${species_name}.all_CDS.fa;
gzip ${OUTDIR}${primate}/${primate}.all_CDS.fa;
