#!/bin/bash
#SBATCH --array=2-19729
#SBATCH --job-name=AlignClust
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
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Clean_clusters/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Alignments_refs/Input_clusters/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_fixed_annotations_refs/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Traits_counts_refs.tsv
mkdir -p ${OUTDIR}

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}SubsetBestOrthos_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${gene_name}
tar xvzf ${INDIR}${gene_name}/${gene_name}.tar.gz --directory ${INDIR}${gene_name}/
rm ${INDIR}${gene_name}/${gene_name}.tar.gz
for filepath in $(ls ${INDIR}${gene_name}/*clean.pep.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1) 
python ${SRC}SubsetBestOrthos.py ${filepath} \
${SPECIES_IDs}summary_species.txt \
${OUTDIR}${gene_name}/${species_name}.clusters.pep.gz
done
#tar arxives again
cd ${INDIR}${gene_name}/
tar cvzf ${gene_name}.tar.gz *clean.pep.gz
rm -r *clean.pep.gz
cd ${OUTDIR}${gene_name}
tar cvzf ${OUTDIR}${gene_name}/${gene_name}.tar.gz *.clusters.pep.gz
rm *.clusters.pep.gz
