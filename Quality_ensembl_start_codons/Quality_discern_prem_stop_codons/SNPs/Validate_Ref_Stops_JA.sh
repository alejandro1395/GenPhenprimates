#!/bin/bash
#SBATCH --array=1-15995
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
TREE_NAMES=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Trees/names.txt
REFERENCE_NAMES=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/summary_species.txt
POSITIONS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_discern_prem_stop_codons/SNPs/Within_stop_codons/
ALIGNMENT=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS/Alignment_CDS/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_discern_prem_stop_codons/SNPs/Between_stop_codons/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_discern_prem_stop_codons/SNPs/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
python ${SRC}Validate_Ref_Stops.py ${ALIGNMENT}${gene_name}/${species_name}.renamed.cds.aln \
${POSITIONS}${gene_name}/All_stop_codons_used_ref.txt \
${gene_name} \
${OUTDIR}${gene_name}/${species_name}.validated.stop_codons_used_ref.txt
