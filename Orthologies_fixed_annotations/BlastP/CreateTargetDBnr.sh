#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/BLAST_in/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/fixed_annotation_results/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_fixed_annotations/BlastP/
mkdir -p ${OUTDIR}BLAST_nrDB
mkdir -p ${OUTDIR}BLAST_nrDB/qu
mkdir -p ${OUTDIR}BLAST_nrDB/out

#Loop for human genes
for filepath in $(ls ${INDIR});
do gene_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}BLAST_nrDB/${gene_name}
#loop for species;
tar tzf ${INDIR}${gene_name}/${gene_name}.tar.gz | while IFS= read -r f ; do 
if [[ $f = *.gz ]] ;
   then echo ">>> Processing file $f"
   tar Oxzf ${INDIR}${gene_name}/${gene_name}.tar.gz "$f" | cat - >> ${OUTDIR}BLAST_nrDB/${gene_name}/allspeciesDB.pep.fa.gz;
fi
done
done
