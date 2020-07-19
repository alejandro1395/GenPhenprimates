#!/bin/bash
#SBATCH --array=1-15
#SBATCH --job-name=QualCDS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=06:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load EMBOSS

#Define PATH argument
REFS_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/RESEQUENCED/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_resequenced/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_resequenced/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1)
python ${SRC}CountStopCodons.py ${INDIR}${gene_name} \
${SPECIES_IDs}pgdp_ids.csv \
${INDIR}${gene_name}/${gene_name}
