#!/bin/bash
#SBATCH --array=1-15671
#SBATCH --job-name=AlignClust
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=00:30:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

#Define PATH argument
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/MUSCLE/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS/SortedSequences/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS/Alignment_proteins/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Genes_predictors/PGLS/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Primate_Traits/OUTPUT/TraitsPerSpecies.txt
mkdir -p ${OUTDIR}

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}SortSeq_repair_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cd ${INDIR}${gene_name}
tar xvzf ${species_name}.tar.gz ${species_name}.renamed.pep.fa.gz
cd .
gunzip -c ${INDIR}${gene_name}/${species_name}.renamed.pep.fa.gz | ${BIN}muscle3.8.31_i86linux64 -in - \
-clwstrictout ${OUTDIR}${gene_name}/${species_name}.renamed.pep.aln
rm ${INDIR}${gene_name}/${species_name}.renamed.pep.*
