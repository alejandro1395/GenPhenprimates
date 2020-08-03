#!/bin/bash
#SBATCH --array=1-15995
#SBATCH --job-name=TransCDS
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
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/MUSCLE/
INDIR_RESEQUENCED=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_resequenced/Check_CDS_frameshift/Resequenced_cds_genes_fasta/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_resequenced/Check_premature_STOP_codon/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_resequenced/Check_premature_STOP_codon/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Primate_Traits/OUTPUT/TraitsPerSpecies.txt
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}Resequenced_pep_genes_fasta/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}Group_by_gene_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}Resequenced_pep_genes_fasta/${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cat ${INDIR_RESEQUENCED}${gene_name}/${species_name}.fasta | \
transeq -sequence stdin -outseq ${OUTDIR}Resequenced_pep_genes_fasta/${gene_name}/${species_name}.qual_raw.pep.fa -trim Y
cat ${OUTDIR}Resequenced_pep_genes_fasta/${gene_name}/${species_name}.qual_raw.pep.fa | sed '/^>/s/.\{2\}$//' > ${OUTDIR}Resequenced_pep_genes_fasta/${gene_name}/${species_name}.qual.pep.fa;
rm ${OUTDIR}Resequenced_pep_genes_fasta/${gene_name}/${species_name}.qual_raw.pep.fa
