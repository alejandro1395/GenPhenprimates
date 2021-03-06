#!/usr/bin/bash
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_in/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_nrDB/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/BlastP/

touch MakeDbBlast_input_JA.txt
for filepath in $(ls ${OUTDIR}*/*.gz);
do echo "$filepath"
echo "gunzip -c ${filepath} | makeblastdb -in - -out ${OUTDIR} -title allspeciesDBnr.pep.fa -dbtype prot" >> MakeDbBlast_input_JA.txt
done
#SUBMISSION TO CLUSTER
#/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${INDIR}out/allspeciesDBnr_100.pep.makeblast.out \
#-e ${INDIR}out/allspeciesDBnr_100.pep.makeblast.err -n makeblastDB -u 1 -t 1 -w 00:05:00
