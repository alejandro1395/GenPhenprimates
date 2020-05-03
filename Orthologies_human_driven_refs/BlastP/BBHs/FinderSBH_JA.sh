#!/bin/bash
#SBATCH --array=1-14
#SBATCH --job-name=FinderSBH
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
INDIR_BLAST=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_out/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/SBHs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/BlastP/BBHs/
mkdir -p ${OUTDIR}

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}FinderSBH_repair_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${gene_name}
tar tzf ${INDIR_BLAST}${gene_name}/${gene_name}.tar.gz | while IFS= read -r f ; do
if [[ $f = *.gz ]] ;
   then echo ">>> Processing file $f"
species_name=$(echo $f | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
   tar Oxzf ${INDIR_BLAST}${gene_name}/${gene_name}.tar.gz "$f" | gunzip -c - | python ${SRC}FinderSBH.py \
${SPECIES_IDs}summary_species.txt \
${species_name} \
${OUTDIR}${gene_name}/${species_name}.pep.SBH.gz
fi
done
#tar arxives again
tar cvzf ${OUTDIR}${gene_name}/${gene_name}.tar.gz ${OUTDIR}${gene_name}/*
rm ${OUTDIR}${gene_name}/*.pep.SBH.gz
