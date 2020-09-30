#!/bin/bash
#SBATCH --array=2-19729
#SBATCH --job-name=OrthoClust
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
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Orthology_clusters/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/Clean_clusters/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_fixed_annotations/BlastP/Clustering/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Traits_counts_refs.tsv
mkdir -p ${OUTDIR}

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}Clean_clusters_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${gene_name}
tar xvzf ${INDIR}${gene_name}/${gene_name}.tar.gz --directory ${INDIR}${gene_name}/
rm ${INDIR}${gene_name}/${gene_name}.tar.gz
for filepath in $(ls ${INDIR}${gene_name}/*.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1) 
python ${SRC}Clean_clusters.py ${filepath} \
${SPECIES_IDs}summary_species.txt \
${OUTDIR}${gene_name}/${species_name}.Clusters_clean.pep.gz
done
#tar arxives again
cd ${INDIR}${gene_name}/
tar cvzf ${gene_name}.tar.gz *.clust.pep.gz
rm *.clust.pep.gz
cd ${OUTDIR}${gene_name}
tar cvzf ${OUTDIR}${gene_name}/${gene_name}.tar.gz *clean.pep.gz
rm *clean.pep.gz
cd .
