#!/bin/bash
#SBATCH --array=1-15
#SBATCH --job-name=PrepCDS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=06:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

#Define PATH argument
REFS_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
REFS_FASTA=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Input_fastas/
RESEQUENCED_FASTA=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_resequenced/MergedAA_resequenced/
RESEQUENCED_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/RESEQUENCED/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_resequenced/TotalPepFasta/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_human_driven_nonrefs/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1)
mkdir -p ${OUTDIR}${gene_name}
for filepath in $(ls ${REFS_FASTA}${gene_name}/*.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cat ${filepath} ${RESEQUENCED_FASTA}${gene_name}/${gene_name}.resequenced.pep.fa.gz > ${OUTDIR}${gene_name}/${species_name}.duplicated.pep.fa.gz
python ${SRC}MergeTotalPep.py ${OUTDIR}${gene_name}/${species_name}.duplicated.pep.fa.gz \
${REFS_IDs}summary_species.txt \
${OUTDIR}${gene_name}/${species_name}.pep.fa.gz \
${RESEQUENCED_IDs}used_reference_species.csv;
rm ${OUTDIR}${gene_name}/${species_name}.duplicated.pep.fa.gz
done
