#!/bin/bash
#SBATCH --array=1-15
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
REFERENCE_GFF=/scratch/devel/lkuderna/PGDP/GENE_VARIANTS_testings/alejandro_testcases/Saguinus_midas.gff
IND_INDELS=/scratch/devel/lkuderna/PGDP/GENE_VARIANTS_testings/alejandro_testcases/PD_0121.variable.filtered.HF.indels.vcf.gz
IND_SNPS=/scratch/devel/lkuderna/PGDP/GENE_VARIANTS_testings/alejandro_testcases/PD_0121.variable.filtered.HF.snps.vcf.gz
SPLITTED_CDS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/Consensus_CDS/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/Concatenated_CDS/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_resequenced/Lukas_tests/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}samples.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
ind_name=$(echo $ARGUMENT1)
mkdir -p ${OUTDIR}${ind_name}
python ${SRC}ConcatenateCDS.py ${SPLITTED_CDS}${ind_name} \
${REFERENCE_GFF} \
chizhang_R002020 \
${OUTDIR}${ind_name}

