#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_in/
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_out/

#We keep the species names for each one of the primates used for annotation from their path
Primates=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/
HUMAN=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/Homo_sapiens/Homo_sapiens.ref_genes.tsv
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_in/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_fixed_annotations/PepAnnotations/

echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}CreateFastaGeneFamilies.py $HUMAN \
${Primates}/summary_species.txt \
$Primates \
$OUTDIR" > ${SRC}qu/CreateFastaGeneFamilies.pep.sh
jobname=$(echo ${SRC}qu/CreateFastaGeneFamilies.pep.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${SRC}out/CreateFastaGeneFamilies.pep.out \
-e ${SRC}out/CreateFastaGeneFamilies.pep.err -n Gene_fam -u 1 -t 1 -w 1-23:59:00 -r lowprio
