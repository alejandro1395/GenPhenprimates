#!/bin/bash
#SBATCH --array=1-2
#SBATCH --job-name=VariantsCalled
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=06:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load VCFTOOLS/0.1.12b
module load TABIX/0.2.5
module load xz/5.2.2
module load BCFTOOLS/1.6

#Define PATH argument
IND_INFO=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/SeparatedVariantsVCFs/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/AllVariantsVCF/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_resequenced/Lukas_tests/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=2 'NR==task_id' ${SRC}TestIndividuals.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
ind_name=$(echo $ARGUMENT1)
bcftools reheader --samples samples.txt -o ${IND_INFO}${ind_name}.indels_renamed.vcf.gz ${IND_INFO}${ind_name}.indels.vcf.gz 
tabix -p vcf ${IND_INFO}${ind_name}.snps.vcf.gz
tabix -p vcf ${IND_INFO}${ind_name}.indels_renamed.vcf.gz
vcf-concat ${IND_INFO}${ind_name}.snps.vcf.gz ${IND_INFO}${ind_name}.indels_renamed.vcf.gz | bgzip -c > ${OUTDIR}${ind_name}.all_var.vcf.gz
vcf-sort ${OUTDIR}${ind_name}.all_var.vcf.gz | bgzip -c > ${OUTDIR}${ind_name}.all_var_sorted.vcf.gz
