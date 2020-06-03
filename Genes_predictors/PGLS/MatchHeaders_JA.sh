#!/bin/bash
#SBATCH --array=1-13
#SBATCH --job-name=MatchHeaders
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
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/PGLS/SortedSequences/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Genes_predictors/PGLS/
TREE=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Trees/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}SortSeq_repair_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
python ${SRC}MatchHeaders.py ${INDIR}${gene_name}/${species_name}.sorted.cds.fa.gz \
${INDIR}${gene_name}/${species_name}.sorted.pep.fa.gz \
${SPECIES_IDs}summary_species.txt \
${TREE}names.txt \
${INDIR}${gene_name}/${species_name}.renamed.cds.fa.gz \
${INDIR}${gene_name}/${species_name}.renamed.pep.fa.gz
#tar file
cd ${INDIR}${gene_name}
tar cvzf ${species_name}.tar.gz ${species_name}.*.fa.gz
cd .
rm ${INDIR}${gene_name}/${species_name}.renamed.*.fa.gz


