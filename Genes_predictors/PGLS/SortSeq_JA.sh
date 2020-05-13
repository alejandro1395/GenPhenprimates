#!/bin/bash
#SBATCH --array=1-14
#SBATCH --job-name=SortCDS
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
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/trimAl/source/trimal
PROT_ALIGNMENT=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Prot_alignments_filtered/
PROT_FASTA=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Input_fastas/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS/SortedSequences/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Genes_predictors/PGLS/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Primate_Traits/OUTPUT/TraitsPerSpecies.txt

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==8' ${SRC}SortSeq_repair_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${gene_name}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cd ${PROT_ALIGNMENT}${gene_name}
tar xvzf ${species_name}.filter2.tar.gz ${species_name}.filter2.90.pep.aln
cd .
python ${SRC}SortSeq.py ${PROT_ALIGNMENT}${gene_name}/${species_name}.filter2.90.pep.aln \
${PROT_FASTA}${gene_name}/${species_name}.pep.fa.gz \
${SPECIES_IDs}summary_species.txt \
${SPECIES_IDs} \
${OUTDIR}${gene_name}/${species_name}.sorted.cds.fa.gz \
${OUTDIR}${gene_name}/${species_name}.sorted.pep.fa.gz
