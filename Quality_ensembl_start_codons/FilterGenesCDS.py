#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
#import seaborn as sns
from scipy import stats
import pickle
from collections import Counter
import copy
from scipy.stats import sem, t
from scipy import mean
import re
import os
import gzip
import time
from Bio import SeqIOHomo_sapiens.ref_genes.tsv
from Bio.Seq import Seq


"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 4:
    human_genes = sys.argv[2]
    reference_gff = sys.argv[1]
    out_file = sys.argv[3]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")



#SPECIES NAMES ASSOCIATES TO TAGS

"""
Read human genes
"""

human_models = pd.read_csv(human_genes, sep='\t', low_memory=False)#panda creation

#DICT SPECIES TAGS
ensembl_ids = human_models['ensembl_peptide_id'].to_list()
gene_names = human_models['hgnc_symbol'].to_list()
ensembl_ids_dict = {ensembl_ids[i]:gene_names[i] for i in range(0, len(ensembl_ids))}


"""
FUNCTIONS
"""


"""
Reverse complement FUNCTIONS
"""

def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)


def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = reverse(seq)
    seq = complement(seq)
    return seq


#MAIN SCRIPT


with gzip.open(reference_gff, "rt") as in_gff:
    last_line = in_gff.readlines()[-1]

total_gene_seq = ""
ref_gene_seq = ""
with gzip.open(reference_gff, "r") as in_gff, \
     open(out_file+"/"+gene_name+"_coordinates.txt", "w") as out_fh:
    for line in in_gff:
        fields = line.rstrip().split("\t")
        line = line.rstrip()
        chrom = fields[0]
        Parent_info = fields[8]
        seq_type = fields[2]
        if seq_type == "CDS" and gene_name in Parent_info:
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            long_seq_record = records[chrom]
            long_seq = long_seq_record.seq
            ref_seq = str(long_seq)[start-1:end]
            ref_gene_seq += ref_seq
            print(chrom+":"+str(start)+"-"+str(end), file=out_fh)
            """
            for position in range(start, end+1):
                new_ident = chrom + "_" + str(position)
                if new_ident in VCF_snps_dict:
                    genotype = get_genotype(VCF_snps_dict[new_ident])
                    total_gene_seq += genotype
                elif new_ident in VCF_indels_dict:
                    genotype = get_genotype(VCF_indels_dict[new_ident])
                    total_gene_seq += genotype
                else:
                    long_seq_record = records[chrom]
                    long_seq = long_seq_record.seq
                    short_seq = str(long_seq)[start-1:end]
                    nt = str(short_seq)[position-start]
                    total_gene_seq += nt

#Obtain final sequence for the gene
    if strand == "-":
        final_seq = reverse_complement(total_gene_seq)
        final_ref_seq = reverse_complement(ref_gene_seq)
    else:
        final_seq = total_gene_seq
        final_ref_seq = ref_gene_seq
    with open(out_file+"/"+gene_name, "w") as out_fh:
        print(">"+gene_name, file=out_fh)
        print(final_seq, file=out_fh)
        print(">ref_"+gene_name, file=out_fh)
        print(final_ref_seq, file=out_fh)
        exit()
        """
