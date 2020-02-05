#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/Query/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastP/
mkdir -p ${OUTDIR}nrDB
mkdir -p ${OUTDIR}nrDB/qu
mkdir -p ${OUTDIR}nrDB/out

touch ${OUTDIR}nrDB/allspeciesDB.pep.gz
for filepath in $(ls ${INDIR}*_*/*.pep.*);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
cat $filepath >> ${OUTDIR}nrDB/allspeciesDB.pep.gz
done

echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python ${SRC}Filternrpep.py ${OUTDIR}nrDB/allspeciesDB.pep.gz ${OUTDIR}nrDB/allspeciesDBnr.pep.gz" > ${OUTDIR}nrDB/qu/allspeciesDBnr.pep.sh
jobname=$(echo ${OUTDIR}nrDB/qu/allspeciesDBnr.pep.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}nrDB/out/allspeciesDBnr.pep.out \
-e ${OUTDIR}nrDB/qu/allspeciesDBnr.pep.err -n merging -u 1 -t 1 -w 02:00:00
#done
