#!/bin/bash
#SBATCH --job-name=QualCDS
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=23:59:00

#Define modules
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load EMBOSS

#Define PATH argument
QUAL_REGIONS=/scratch/devel/lkuderna/PGDP/project_NIST/benchmark/toplevel_id_refmt.MASK.bed
HUMAN_GFF=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/Homo_sapiens/Homo_sapiens.gff.gz
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Lukas_tests/CheckHighQualRegions/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Lukas_tests/CheckHighQualRegions/

# EXECUTING PART
python -u ${SRC}CheckRegionsGFF.py ${SRC}human_genes_clusters \
${HUMAN_GFF} \
${QUAL_REGIONS} \
${OUTDIR}non-covered_intervals.txt
