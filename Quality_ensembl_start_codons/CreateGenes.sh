#!/bin/bash
#SBATCH --array=1-15995
#SBATCH --job-name=PrepCDS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=01:00:00

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
Primates="Aotus_nancymaae"
REFS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Fastas/REFS/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/GFFs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_ensembl_start_codons/
HUMAN_MODEL=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/Human_models/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/Ensembl_genes/Aotus_nancymaae/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
python ${SRC}CreateGenes.py ${INDIR}Aotus_nancymaae.cds.fa.gz \
${gene_name} \
${INDIR}Aotus_nancymaae.orthos.txt \
${HUMAN_MODEL}Homo_sapiens.ref_genes.tsv \
${OUTDIR}${gene_name}/${gene_name}.fa \
