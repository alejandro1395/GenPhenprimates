#!/bin/bash
#SBATCH --array=1-16230
#SBATCH --job-name=AlignTrim
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=05:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

#Define PATH argument
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/trimAl/source/trimal
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/Prot_alignments/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/Prot_alignments_filtered/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_fixed_annotations_refs/Trimming/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Traits_counts_refs.tsv
mkdir -p ${OUTDIR}

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}trimAl_filter1_pos_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
dir_out=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${dir_out}
species_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
gunzip ${ARGUMENT1}
input_file=$(echo $ARGUMENT1 | rev | cut -c 4- | rev)
${BIN} -in $input_file \
-gt 0.9 -out ${OUTDIR}${dir_out}/${species_name}.filter1.pep.aln -colnumbering > ${OUTDIR}${dir_out}/${species_name}.filter1.ref_pos.txt
gzip $input_file
cd ${OUTDIR}${dir_out}
tar cvzf ${dir_out}.tar.gz *
rm *filter1*
cd .
