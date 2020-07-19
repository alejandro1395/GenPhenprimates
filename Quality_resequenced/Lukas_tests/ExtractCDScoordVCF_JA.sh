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
REFERENCE_FASTA=/scratch/devel/lkuderna/PGDP/GENE_VARIANTS_testings/alejandro_testcases/Saguinus_midas.fasta
REFERENCE_GFF=/scratch/devel/lkuderna/PGDP/GENE_VARIANTS_testings/alejandro_testcases/Saguinus_midas.gff
IND_INDELS=/scratch/devel/lkuderna/PGDP/GENE_VARIANTS_testings/alejandro_testcases/PD_0121.variable.filtered.HF.indels.vcf.gz
IND_SNPS=/scratch/devel/lkuderna/PGDP/GENE_VARIANTS_testings/alejandro_testcases/PD_0121.variable.filtered.HF.snps.vcf.gz
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/Coordinates_gff/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_resequenced/Lukas_tests/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=11 'NR==task_id' ${SRC}TestGenes.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
gene_name=$(echo $ARGUMENT1)
mkdir -p ${OUTDIR}${gene_name}
python ${SRC}ExtractCDScoordVCF.py ${REFERENCE_FASTA} \
${REFERENCE_GFF} \
${IND_INDELS} \
${IND_SNPS} \
$gene_name \
${OUTDIR}${gene_name}

