#!/bin/bash
#SBATCH --array=1-15995
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
REFERENCE_NAMES=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/summary_species.txt
TREE_NAMES=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Trees/names.txt
ALIGNMENT_PATH=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Prot_alignments_filtered/
INDIVIDUALS_PATH=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_genome_bottle/Individuals_cds_fastas/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_genome_bottle/Individuals_cds_genes_fastas/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_genome_bottle/CDS_sequences/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}Group_by_gene_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cd ${ALIGNMENT_PATH}${gene_name}
tar xvzf ${species_name}.filter2.tar.gz ${species_name}.filter2.90.pep.aln
cd .
python ${SRC}Group_by_gene.py ${INDIVIDUALS_PATH} \
${ALIGNMENT_PATH}${gene_name}/${species_name}.filter2.90.pep.aln \
${REFERENCE_NAMES} \
${TREE_NAMES} \
${OUTDIR}${gene_name}/${species_name}.fasta
