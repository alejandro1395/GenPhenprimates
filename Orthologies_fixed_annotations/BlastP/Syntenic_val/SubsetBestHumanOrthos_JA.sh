#!/bin/bash
#SBATCH --array=1-19728
#SBATCH --job-name=OrthoClust
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
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Orthology_clusters/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/BlastP/Syntenic_val/
TRAITS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Phenomes/Primate_Traits/OUTPUT/TraitsPerSpecies.txt
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Orthology_synteny_val/
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}qu
mkdir -p ${OUTDIR}out

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}Clean_clusters_repair_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${gene_name}
tar xvzf ${INDIR}${gene_name}/${gene_name}.tar.gz --directory ${INDIR}${gene_name}/
for filepath in $(ls ${INDIR}${gene_name}${INDIR}${gene_name}/*clean.pep.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1) 
cd ${INDIR}${gene_name}
tar xvzf ${gene_name}.tar.gz ${INDIR}${gene_name}/${species_name}.Clusters_clean.pep.gz
cd .
python ${SRC}SubsetBestHumanOrthos.py ${INDIR}${gene_name}/${species_name}.Clusters_clean.pep.gz \
${SPECIES_IDs}summary_species.txt \
${OUTDIR}${gene_name}/${species_name}.Validation_clusters.pep.gz
rm ${INDIR}${gene_name}/${species_name}.Clusters_clean.pep.gz
done
