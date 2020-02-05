#!/bin/bash

module load purge
module load PYTHON/3.6.3
module load BLAST+

mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_int
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_out

#We keep the species names for each one of the primates used for annotation from their path
Primates=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_int



DATA=/scratch/devel/avalenzu/Impute_Master_Project/ANALYSIS_jan2019_panel300/BEAGLE_4.0_analysis_chimp/SAMPLES_vcf/
MAP=/scratch/devel/avalenzu/Impute_Master_Project/ANALYSIS_jan2019_panel300/BEAGLE_4.1_analysis_chimp/MAP/
BIN=/scratch/devel/avalenzu/Impute_Master_Project/ANALYSIS_jan2019_panel300/bin/
REF=/scratch/devel/avalenzu/Impute_Master_Project/ANALYSIS_jan2019_panel300/BEAGLE_4.0_analysis_chimp/PANEL/SPLITTED_VCF/
OUTDIR=/scratch/devel/avalenzu/Impute_Master_Project/ANALYSIS_jan2019_panel300/BEAGLE_4.0_analysis_chimp/IMPUTATION/

mkdir -p ${OUTDIR}/qu
mkdir -p ${OUTDIR}/out
count=0
for filepath in $(ls ${REF}filtered*.vcf.gz);
do count=$(($count+1))
echo $count

echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load xz/5.2.2
module load SAMTOOLS/1.3
module load java
module load GATK/4.0.8.1
module load TABIX/0.2.6
module load VCFTOOLS/0.1.7

java -Xmx70g -Djava.io.tmpdir=${OUTDIR}tmp/ \
-XX:-UseGCOverheadLimit \
-jar ${BIN}beagle.r1399.jar \
gl=${DATA}VCF_Boe1_merged_nonmiss.vcf.recode.vcf \
ref=${REF}filtered_chimp_chr21_Panel_splitted.${count}.vcf.gz \
nthreads=8 \
out=${OUTDIR}VCF_Boe1_merged_nonmiss_${count} \
map=${MAP}final_chr21.map \
impute=true \
gprobs=true" > ${OUTDIR}qu/VCF_Boe1_merged_gl_imp_${count}.sh
jobname=$(echo ${OUTDIR}qu/VCF_Boe1_merged_gl_imp_${count}.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/VCF_Boe1_merged_gl_imp_${count}.out \
-e ${OUTDIR}out/VCF_Boe1_merged_gl_imp_${count}.sh -n gl_imp_${count} -u 8 -t 1 -w 20:00:00


#gl=../../SAMPLES_LOW/VCF_Boe1_merged/VCF_Boe1_merged_nonmiss.vcf.gz
done
