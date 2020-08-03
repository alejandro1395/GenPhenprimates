#!/bin/bash
#SBATCH --array=1-343
#SBATCH --job-name=PrepIndCDS
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
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_resequenced/Check_CDS_frameshift/
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}Individuals_cds_fastas/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_resequenced/Check_CDS_frameshift/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}Individuals_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
individual_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
python ${SRC}Individuals_to_fasta.py $ARGUMENT1 \
${OUTDIR}Individuals_cds_fastas/${individual_name}.fasta
