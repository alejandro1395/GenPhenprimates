#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_in/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/BLAST_nrDB
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Orthologies_human_driven_refs/BlastP/
mkdir -p ${OUTDIR}BLAST_nrDB
mkdir -p ${OUTDIR}BLAST_nrDB/qu
mkdir -p ${OUTDIR}BLAST_nrDB/out

#Loop for human genes
for filepath in $(ls ${INDIR});
do gene_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
echo $gene_name
mkdir -p ${OUTDIR}BLAST_nrDB/${gene_name}
#loop for species;
tar tzf ${INDIR}${gene_name}/${gene_name}.tar.gz | while IFS= read -r f ; do
if [[ $f = *.gz ]] ;
   then echo ">>> Processing file $f"
   tar Oxzf ${INDIR}${gene_name}/${gene_name}.tar.gz "$f" | gunzip -c - | blastp -query - \
-db ${TARGETDB}allspeciesDBnr_100.pep \
-evalue 1e-4 \
-num_alignments 0 \
-out ${OUTDIR}${species_name}/${species_name}.nr100.pep \
-num_threads 4;
gzip ${OUTDIR}${species_name}/${species_name}.nr100.pep
 >> ${OUTDIR}BLAST_nrDB/${gene_name}/allspeciesDB.pep.fa.gz;
fi
done
done



for filepath in $(ls ${QUERY}*/*.nr100.pep.fa.gz);
do species_name=$(echo $filepath | rev | cut -d'/' -f1 | rev | cut -d \. -f 1)
mkdir -p ${OUTDIR}${species_name}/
echo "#!/bin/bash
module purge
module unload gcc/4.9.3-gold
module load gcc/6.3.0
module load PYTHON/3.6.3
module load BLAST+

gunzip -c ${filepath} | blastp -query - \
-db ${TARGETDB}allspeciesDBnr_100.pep \
-evalue 1e-4 \
-num_alignments 0 \
-out ${OUTDIR}${species_name}/${species_name}.nr100.pep \
-num_threads 4;
gzip ${OUTDIR}${species_name}/${species_name}.nr100.pep" > ${OUTDIR}qu/${species_name}.nr100.pep.sh
jobname=$(echo ${OUTDIR}qu/${species_name}.nr100.pep.sh)
chmod 755 $jobname

#SUBMISSION TO CLUSTER
/scratch/devel/avalenzu/CNAG_interface/submit.py -c ${jobname} -o ${OUTDIR}out/${species_name}.nr100.pep.out \
-e ${OUTDIR}out/${species_name}.nr100.pep.err -n ${species_name} -u 4 -t 1 -w 1-23:55:00 -r lowprio
done
