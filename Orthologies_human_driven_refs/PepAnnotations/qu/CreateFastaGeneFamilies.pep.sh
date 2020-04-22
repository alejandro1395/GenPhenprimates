#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/PepAnnotations/CreateFastaGeneFamilies.py /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/Homo_sapiens/Homo_sapiens.ref_genes.tsv /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS//summary_species.txt /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/ /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_in/
