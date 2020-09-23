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
FASTA_CDS_RESEQUENCED=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_genome_bottle/Individuals_cds_genes_fastas/
FASTA_CDS_REFS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS/SortedSequences/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_genome_bottle/Alignments_within_ref/
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}Fastas/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_genome_bottle/Alignments_within_ref/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}Fastas/${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cd ${FASTA_CDS_REFS}${gene_name}
tar xvzf ${species_name}.tar.gz ${species_name}.renamed.cds.fa.gz
cd .
python ${SRC}PrepareCDS.py ${FASTA_CDS_REFS}${gene_name}/${species_name}.renamed.cds.fa.gz \
${REFERENCE_NAMES} \
${TREE_NAMES} \
${FASTA_CDS_RESEQUENCED}${gene_name}/${species_name}.fasta \
${OUTDIR}Fastas/${gene_name}/ \
${gene_name}
rm ${FASTA_CDS_REFS}${gene_name}/${species_name}.renamed.cds.fa.gz
