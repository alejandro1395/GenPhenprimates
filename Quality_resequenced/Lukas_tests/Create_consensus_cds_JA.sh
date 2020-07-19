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
REFERENCE_FASTA=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/Ref_fasta/Saguinus_midas.fasta
IND_VCF=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/AllVariantsVCF/
IND_COORD=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/Coordinates_gff/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/Consensus_CDS/
mkdir -p ${OUTDIR}
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_resequenced/Lukas_tests/

# Define arguments in each task
ARGUMENT1=`awk -v task_id=1 'NR==task_id' ${SRC}samples.txt`
# Print info of the task
echo $ARGUMENT1

# EXECUTING PART
ind_name=$(echo $ARGUMENT1)
mkdir -p ${OUTDIR}${ind_name}
count=0
cat ${IND_COORD}${ind_name}/*_coordinates.txt | while read line;
do count=$((count + 1)) 
samtools faidx ${REFERENCE_FASTA} ${line} | bcftools consensus --haplotype R ${IND_VCF}${ind_name}.variable.filtered.HF.all_var_sorted.vcf.gz > ${OUTDIR}${ind_name}/${ind_name}.cds.chunk${count}.fa
done
