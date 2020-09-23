#!/bin/bash
#SBATCH --array=1-15995
#SBATCH --job-name=PrepCDS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=01:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

#Define PATH argument
TREE_NAMES=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Trees/names.txt
REFERENCE_NAMES=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/summary_species.txt
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_genome_bottle/Alignments_within_ref/Alignments/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_genome_bottle/SNPs/
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}Within_stop_codons/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_genome_bottle/SNPs/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}Within_stop_codons/${gene_name}
tar xvzf ${INDIR}${gene_name}/${gene_name}.tar.gz --directory ${INDIR}${gene_name}/
for filepath in $(ls ${INDIR}${gene_name}/*within.cds.aln);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
python ${SRC}Count_Ref_Stops.py ${INDIR}${gene_name}/${species_name}.within.cds.aln \
${gene_name} \
${species_name} \
${OUTDIR}Within_stop_codons/${gene_name}/${species_name}.stop_codons.txt
done
#tar arxives again
rm ${INDIR}${gene_name}/${gene_name}.tar.gz
cd ${INDIR}${gene_name}/
tar cvzf ${gene_name}.tar.gz *.within.cds.aln
cd .
rm -r ${INDIR}${gene_name}/*.within.cds.aln
ls ${OUTDIR}Within_stop_codons/${gene_name}/*.stop_codons.txt | while read file; do cat $file | tail -n+2 >> ${OUTDIR}Within_stop_codons/${gene_name}/All_stop_codons.txt; done
rm ${OUTDIR}Within_stop_codons/${gene_name}/*.stop_codons.txt
