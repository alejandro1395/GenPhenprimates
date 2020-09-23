#!/usr/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#Create folders
mkdir -p /scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/BBHs/

#We keep the species names for each one of the primates used for annotation from their path
SPECIES_IDs=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Annotations/REFS/
CHROM_path=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/data/Genomes/Syntenic_alignments/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/results/Orthology_synteny_val/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_refs/BlastP/Syntenic_val/
mkdir -p ${OUTDIR}qu
mkdir -p ${OUTDIR}out

validating_species="Chlorocebus_sabaeus Pan_troglodytes Otolemur_garnettii Pan_paniscus Papio_anubis Pongo_abelii"

for filepath in $(ls ${OUTDIR}*/*Validation_clusters.pep.gz);
do species_name_ref=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
echo $validating_species | tr " " "\n" | while read species2;
do echo "#!/bin/bash
module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3
module load UCSCTOOLS/389

liftOver ${OUTDIR}${species_name_ref}/${species_name_ref}.ref.Homo_sapiens.Clusters_loc.pep.gz \
${CHROM_path}${species2}/*.rbest.chain.gz \
${OUTDIR}${species_name_ref}/${species_name_ref}.ref.${species2}.Lifted.tsv \
${OUTDIR}${species_name_ref}/unMapped_${species2}.tsv" > ${OUTDIR}qu/${species_name_ref}.ref.${species2}.Lifted.pep.sh
jobname=$(echo ${OUTDIR}qu/${species_name_ref}.ref.${species2}.Lifted.pep.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name_ref}.ref.${species2}.Lifted.pep.out \
-e ${OUTDIR}out/${species_name_ref}.ref.${species2}.Lifted.pep.err -n Lift${species_name} -u 1 -t 1 -w 06:00:00
done
done
