#!/bin/bash
#SBATCH --array=1-19729
#SBATCH --job-name=BlastPAll
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

QUERY=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_in/
TARGET=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_nrDB/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_out/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/BlastP/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}BlastPallqueryVSallnr_100_repair_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1 | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${gene_name}
tar tzf ${QUERY}${gene_name}/${gene_name}.tar.gz | while IFS= read -r f ; do
if [[ $f = *.gz ]] ;
   then echo ">>> Processing file $f"
species_name=$(echo $f | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
   tar Oxzf ${QUERY}${gene_name}/${gene_name}.tar.gz "$f" | zcat - | blastp -query - \
-db ${TARGET}${gene_name}/allspeciesDBnr_100.pep \
-evalue 1e-4 \
-num_alignments 0 \
-out ${OUTDIR}${gene_name}/${species_name}.nr100.pep \
-num_threads 4;
gzip ${OUTDIR}${gene_name}/${species_name}.nr100.pep
fi
done
#tar arxives again
tar cvzf ${OUTDIR}${gene_name}/${gene_name}.tar.gz ${OUTDIR}${gene_name}/*
rm ${OUTDIR}${gene_name}/*.pep.gz
