#!/bin/bash
#SBATCH --array=6-15
#SBATCH --job-name=AlignClust
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=01:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

#Define PATH argument
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/MUSCLE/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS_resequenced/SortedSequences/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS_resequenced/Alignment_proteins/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Genes_predictors_resequenced/PGLS/
mkdir -p ${OUTDIR}

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}SortCDS_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
gunzip -c ${INDIR}${gene_name}/${species_name}.renamed.pep.fa.gz | ${BIN}muscle3.8.31_i86linux64 -in - \
-clwstrictout ${OUTDIR}${gene_name}/${species_name}.renamed.pep.aln
