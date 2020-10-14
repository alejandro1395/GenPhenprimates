#!/bin/bash
#SBATCH --array=1-16202
#SBATCH --job-name=AlignQual
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=05:00:00

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/trimAl/source/trimal
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/Prot_alignments_filtered/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/Prot_alignments_qual/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_fixed_annotations_refs/Qualities/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Traits_counts_refs.tsv

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}../Trimming/trimAl_filter2_spc_input_JA.txt`
# Print info of the task

# EXECUTING PART
dir_out=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${dir_out}
tar xvzf ${INDIR}${dir_out}/${dir_out}.tar.gz --directory ${INDIR}${dir_out}/
rm ${INDIR}${dir_out}/${dir_out}.tar.gz
for filepath in $(ls ${INDIR}${dir_out}/*filter2.90.pep.aln);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
${BIN} -in ${filepath} \
-sident > ${OUTDIR}${dir_out}/${species_name}.SpeciesIdentities
done
#rm tar files
cd ${INDIR}${dir_out}/
tar cvzf ${dir_out}.tar.gz *
rm *.pep.aln
rm *.txt
