#!/usr/bin/bash

#SBATCH -J trimAl_spc_info # A single job name for the array
#SBATCH -u 1
#SBATCH -t 1 # one core
#SBATCH -n 1 # on one node
#SBATCH -t 0-2:00 # Running time of 2 hours
#SBATCH –mem 4000 # Memory request of 4 GB
#SBATCH -o tophat_%A_%a.out # Standard output
#SBATCH -e tophat_%A_%a.err # Standard error
# grab out filename from the array exported from our ‘parent’ shell
FILENAME=${FILES[$SLURM_ARRAY_TASK_ID]}
# make & move into new directory, and run!
mkdir ${FILENAME}_out
cd ${FILENAME}_out
tophat /n/informatics/databases/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex \
../$FILENAME 


/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.${gene_name}.trimAl_spc_info.out \
#-e ${OUTDIR}out/${species_name}.${gene_name}.trimAl_spc_info.err -n ${species_name}.${gene_name} -u 1 -t 1 -w 00:05:00

