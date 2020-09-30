#!/bin/bash

module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_in/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_nrDB/
mkdir -p $OUTDIR
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_fixed_annotations/BlastP/

for filepath in $(ls ${OUTDIR}*/*.gz);
do dir_out=$(echo $(dirname $filepath)) 
gunzip -c $filepath | makeblastdb -in - \
-out ${dir_out}/allspeciesDBnr_100.pep \
-title allspeciesDBnr_100.pep \
-dbtype prot
done
#jobname=$(echo ${INDIR}qu/allspeciesDBnr_100.pep.makeblast.sh)
#chmod 755 $jobname

#SUBMISSION TO CLUSTER
#/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${INDIR}out/allspeciesDBnr_100.pep.makeblast.out \
#-e ${INDIR}out/allspeciesDBnr_100.pep.makeblast.err -n makeblastDB -u 1 -t 1 -w 00:05:00
