#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

python /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_fixed_annotations/PepAnnotations/CreateFastaGeneFamilies.py /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/Homo_sapiens/Homo_sapiens.ref_genes.tsv /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED//summary_species.txt /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS_FIXED/ /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_in/
