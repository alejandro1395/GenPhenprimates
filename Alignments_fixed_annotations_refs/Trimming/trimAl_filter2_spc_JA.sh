#!/bin/bash
#SBATCH --array=1-16202
#SBATCH --job-name=AlignTrimm
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
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/Prot_alignments_filtered/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_fixed_annotations_refs/Trimming/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Traits_counts_refs.tsv

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}trimAl_filter2_spc_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
dir_out=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${dir_out}
tar xvzf ${OUTDIR}${dir_out}/${dir_out}.tar.gz --directory ${OUTDIR}${dir_out}/
rm ${OUTDIR}${dir_out}/${dir_out}.tar.gz
for i in {50..100..2}
do for filepath in $(ls ${OUTDIR}${dir_out}/*filter1.pep.aln);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
${BIN} -in $filepath \
-resoverlap 0.90 -seqoverlap $i > ${OUTDIR}${dir_out}/${species_name}.filter2.${i}.pep.aln
done
done
#rm tar files
cd ${OUTDIR}${dir_out}/
tar cvzf ${dir_out}.tar.gz *
rm *.pep.aln
rm *.txt
cd .
