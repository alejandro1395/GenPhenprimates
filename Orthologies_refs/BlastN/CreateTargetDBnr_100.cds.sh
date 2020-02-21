#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/Query/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastN/
mkdir -p ${OUTDIR}nrDB
mkdir -p ${OUTDIR}nrDB/qu
mkdir -p ${OUTDIR}nrDB/out

#loop for concatenating all reference species information
if [ -f ${OUTDIR}nrDB/allspeciesDB.cds.fa.gz ] ; then
    rm ${OUTDIR}nrDB/allspeciesDB.cds.fa.gz
fi
touch ${OUTDIR}nrDB/allspeciesDB.cds.fa.gz
for filepath in $(ls ${INDIR}*_*/*.nr100.cds.fa.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cat $filepath >> ${OUTDIR}nrDB/allspeciesDB.cds.fa.gz
done

echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3
module load CDHIT/4.7

gunzip ${OUTDIR}nrDB/allspeciesDB.cds.fa.gz;
cd-hit \
-i ${OUTDIR}nrDB/allspeciesDB.cds.fa \
-o ${OUTDIR}nrDB/allspeciesDBnr_100.cds.fa \
-c 1.0 \
-n 5 \
-M 16000 \
-d 0 \
-T 8;
gzip ${OUTDIR}nrDB/allspeciesDB.cds.fa;
if [ -f ${OUTDIR}nrDB/allspeciesDB.cds.fa ] ; then
    rm ${OUTDIR}nrDB/allspeciesDB.cds.fa
fi
gzip ${OUTDIR}nrDB/allspeciesDBnr_100.cds.fa" > ${OUTDIR}nrDB/qu/allspeciesDBnr_100.cds.fa.sh
jobname=$(echo ${OUTDIR}nrDB/qu/allspeciesDBnr_100.cds.fa.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}nrDB/out/allspeciesDBnr_100.cds.out \
-e ${OUTDIR}nrDB/qu/allspeciesDBnr_100.cds.err -n merging -u 1 -t 1 -w 02:O0:00
#done
