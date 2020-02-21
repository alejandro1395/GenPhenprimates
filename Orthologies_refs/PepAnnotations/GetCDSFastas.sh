#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3
module load BEDTools

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_out/

#We keep the species names for each one of the primates used for annotation from their path
Primates="Carlito_syrichta Propithecus_coquereli"
REFS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Fastas/REFS/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/PepAnnotations/
mkdir -p ${SRC}qu/
mkdir -p ${SRC}out/

echo $Primates | tr " " "\n" | while read primate;
#CREATE GFF FILE WITH JUST CDS ANNOTATION 
do zcat ${OUTDIR}${primate}/*.gff.gz | awk -v OFS='\t' -v FS='\t' '$3 == "CDS" {print $0}' > \
${OUTDIR}${primate}/${primate}.CDS.gff;
gzip ${OUTDIR}${primate}/${primate}.CDS.gff;
#IMPORT ZIPPED FILE INTO BEDTOOLS GETFASTA
zcat ${OUTDIR}${primate}/${primate}.CDS.gff.gz | bedtools getfasta -fi ${REFS}${primate}/*.fna -bed - -s -fo stdout | \
paste - - > ${OUTDIR}${primate}/${primate}.all_CDS.fa;
gzip ${OUTDIR}${primate}/${primate}.all_CDS.fa;

#SUBMISSION TO CLUSTER
#/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}Query/out/${species_name}.nr99.pep.out \
#-e ${OUTDIR}Query/out/${species_name}.nr99.pep.err -n ${species_name} -u 1 -t 1 -w 00:05:00
done
