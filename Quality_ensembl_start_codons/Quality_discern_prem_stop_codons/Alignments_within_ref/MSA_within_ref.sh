#!/bin/bash
#SBATCH --array=2-15995
#SBATCH --job-name=AlignClust
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=02:00:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+
module load EMBOSS

#Define PATH argument
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
BIN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/bin/MUSCLE/
ALIGNMENT=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_discern_prem_stop_codons/Alignments_within_ref/Alignments/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/Ensembl_genes/Aotus_nancymaae/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/Alignments/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_ensembl_start_codons/Quality_discern_prem_stop_codons/Alignments_within_ref/
mkdir -p ${OUTDIR}

# Define arguments in each task
ARGUMENT1=`awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}PrepareCDS_input_JA.txt`
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo $ARGUMENT1


# EXECUTING PART
species_name="Aotus_nancymaae"
gene_name=$(echo $(basename $(dirname $ARGUMENT1)))
mkdir -p ${OUTDIR}${gene_name}
tar xvzf ${ALIGNMENT}${gene_name}/${gene_name}.tar.gz --directory ${ALIGNMENT}${gene_name}/
cat ${ALIGNMENT}${gene_name}/${species_name}.within.cds.aln | \
seqret -sequence stdin -osformat2 fasta -outseq stdout | \
${BIN}muscle3.8.31_i86linux64 -profile -in1 - \
-in2 ${INDIR}${gene_name}/${gene_name}.fa -clwstrictout ${OUTDIR}${gene_name}/${species_name}.within.cds.aln
#tar arxives again
rm ${ALIGNMENT}${gene_name}/${gene_name}.tar.gz
cd ${ALIGNMENT}${gene_name}/
tar cvzf ${gene_name}.tar.gz *.within.cds.aln
cd .
rm -r ${ALIGNMENT}${gene_name}/*.within.cds.aln
