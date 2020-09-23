#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
Primates="Aotus_nancymaae"
REFS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Fastas/REFS/
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/GFFs/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Quality_ensembl_start_codons/
HUMAN_MODEL=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Quality_ensembl_start_codons/Human_models/

echo $Primates | tr " " "\n" | while read primate;
do python ${SRC}MergeCDSFastas.py ${INDIR}${primate}.all_CDS.fa.gz ${INDIR}${primate}.CDS.gff.gz \
${INDIR}Aotus_nancymaae.orthos.txt ${HUMAN_MODEL}Homo_sapiens.ref_genes.tsv ${INDIR}${primate}.cds.fa.gz;
#SUBMISSION TO CLUSTER
#/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.nr.pep.out \
#-e ${OUTDIR}out/${species_name}.nr.pep.err -n ${species_name} -u 1 -t 1 -w 00:05:00
done
