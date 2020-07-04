#!/bin/bash
#SBATCH --array=1-772
#SBATCH --job-name=Cons
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=05:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load EMBOSS

#Define PATH argument
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/RESEQUENCED/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/MUSCLE/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Consensus_alignments/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Consensus_seq/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Consensus_resequenced/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}Consensus_seq_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
dir_out=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${dir_out}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cons -sequence $ARGUMENT1 -outseq ${OUTDIR}${dir_out}/${species_name}.fa -name ${species_name};
gzip ${OUTDIR}${dir_out}/${species_name}.fa
