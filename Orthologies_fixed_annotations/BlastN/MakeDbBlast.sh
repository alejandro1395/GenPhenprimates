#!/bin/bash

#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BLAST_in/nrDB/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastN/

echo "#!/bin/bash
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

gunzip -c ${INDIR}allspeciesDBnr_100.cds.fa.gz | makeblastdb -in - \
-out ${INDIR}allspeciesDBnr_100.cds \
-title allspeciesDBnr_100.cds.fa \
-dbtype nucl" > ${INDIR}qu/allspeciesDBnr_100.cds.makeblast.sh
jobname=$(echo ${INDIR}qu/allspeciesDBnr_100.cds.makeblast.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${INDIR}out/allspeciesDBnr_100.cds.makeblast.out \
-e ${INDIR}out/allspeciesDBnr_100.cds.makeblast.err -n makeblastDB -u 1 -t 1 -w 00:05:00
