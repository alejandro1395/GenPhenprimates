#!/bin/bash
#SBATCH --array=6-15
#SBATCH --job-name=SortCDS
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
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Trees/
CDS_REF=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
CDS_RESEQUENCED=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_resequenced/MergedCDS_resequenced/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/trimAl/source/trimal
PROT_ALIGNMENT=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_resequenced/Prot_alignments_filtered/
PROT_REF=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Input_fastas/
PROT_RESEQUENCED=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_resequenced/TotalPepFasta/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS_resequenced/SortedSequences/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Genes_predictors_resequenced/PGLS/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}SortCDS_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
python ${SRC}SortSeq.py ${PROT_ALIGNMENT}${gene_name}/${species_name}.filter2.90.pep.aln \
${PROT_REF}${gene_name}/${species_name}.pep.fa.gz \
${PROT_RESEQUENCED}${gene_name}/${species_name}.pep.fa.gz \
${SPECIES_IDs}names.txt \
${CDS_REF} \
${CDS_RESEQUENCED}${gene_name}/${gene_name}.resequenced.cds.fa.gz \
${OUTDIR}${gene_name}/${species_name}.sorted.cds.fa.gz \
${OUTDIR}${gene_name}/${species_name}.sorted.pep.fa.gz
