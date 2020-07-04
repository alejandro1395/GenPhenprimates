#!/bin/bash
#SBATCH --array=1-15
#SBATCH --job-name=PrepCDS
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
REFS_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/RESEQUENCED/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Consensus_resequenced/
CONSENSUS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Consensus_seq/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_resequenced/MergedCDS_resequenced/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_human_driven_nonrefs/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1)
mkdir -p ${OUTDIR}${gene_name}
for filepath in $(ls ${INDIR}${gene_name}/*.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
python ${SRC}PrepareCDS_noncons.py ${filepath} \
${gene_name} \
${species_name} \
${CONSENSUS}${gene_name} \
${OUTDIR}${gene_name}/ \
${SPECIES_IDs}/to_exclude_ids.csv
done
