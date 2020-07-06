#!/bin/bash
#SBATCH --array=2-15
#SBATCH --job-name=PrepCDS
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
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_resequenced/MergedCDS_resequenced/
CONSENSUS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Consensus_seq/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_resequenced/MergedAA_resequenced/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_human_driven_nonrefs/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1)
mkdir -p ${OUTDIR}${gene_name}
for filepath in $(ls ${INDIR}${gene_name}/*.gz);
do zcat ${filepath} | \
transeq -sequence stdin -sformat pearson -outseq ${OUTDIR}${gene_name}/${gene_name}.resequenced_raw.pep.fa -trim Y -clean Y
cat ${OUTDIR}${gene_name}/${gene_name}.resequenced_raw.pep.fa | sed '/^>/s/.\{2\}$//' > ${OUTDIR}${gene_name}/${gene_name}.resequenced.pep.fa;
rm ${OUTDIR}${gene_name}/${gene_name}.resequenced_raw.pep.fa
gzip ${OUTDIR}${gene_name}/${gene_name}.resequenced.pep.fa
done
