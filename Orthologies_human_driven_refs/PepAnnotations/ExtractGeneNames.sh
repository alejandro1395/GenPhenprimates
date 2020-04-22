#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastP/

#loop for concatenating all reference species information
for filepath in $(ls ${INDIR}*/*.gff.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
zcat $filepath | grep genename | cut -f 9 | tr "=" "\t" | tr ";" "\t" | cut -f 2,4 > ${INDIR}${species_name}/${species_name}.ref_genes.tsv
echo $filepath processed;
done

#FOR PROPITHECUS AND CARLITO
#zcat Carlito_syrichta.CDS.gff.gz | grep gene | tr ";" "\t" | awk '{printf "%s",$9; for (i=1;i<=NF;i++) { if ($i ~ /gene/) { printf " %s",$i; }  } print ""  }'
#cat Carlito_syrichta.ref_genes_raw.tsv | tr -d "ID=" | tr -d "gene=" | tr " " "\t" | sort | uniq
