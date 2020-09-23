#!/bin/bash
#SBATCH --array=1-2
#SBATCH --job-name=MakeDBBlast
#SBATCH --output=out/MakeDBBlast_%A_%a.out
#SBATCH --error=qu/MakeDBBlast_%A_%a.err
#SBATCH --time=01:00:00

#Define PATH argument

SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/BlastP/

# Define arguments in each task
arguments='awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id' ${SRC}MakeDbBlast_input_JA.txt'
echo $SLURM_ARRAY_TASK_ID
# Print info of the task
echo ${arguments}

# Execute
#gunzip -c ${SRC}MakeDbBlast.sh ${arguments}
