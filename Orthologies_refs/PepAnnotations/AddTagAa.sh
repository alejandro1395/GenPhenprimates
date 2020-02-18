#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3
module load EMBOSS

#We keep the species names for each one of the primates used for annotation from their path
Primates="Carlito_syrichta Propithecus_coquereli"
REFS=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Fastas/REFS/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/PepAnnotations/

echo $Primates | tr " " "\n" | while read primate;
do zcat ${OUTDIR}${primate}/${primate}.pep.fa.gz | perl -p -e 's/^>(.*)/>'$primate' $1/' > ${OUTDIR}${primate}/${primate}.pep.fa;
rm ${OUTDIR}${primate}/${primate}.pep.fa.gz;
gzip ${OUTDIR}${primate}/${primate}.pep.fa;

#SUBMISSION TO CLUSTER
#/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.nr.pep.out \
#-e ${OUTDIR}out/${species_name}.nr.pep.err -n ${species_name} -u 1 -t 1 -w 00:05:00
done
