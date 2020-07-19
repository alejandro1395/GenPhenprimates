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
from Bio import SeqIO
from Bio.Seq import Seq


"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 7:
    reference_fasta = sys.argv[1]
    reference_gff = sys.argv[2]
    ind_indels = sys.argv[3]
    ind_snps = sys.argv[4]
    gene_name = sys.argv[5]
    out_file = sys.argv[6]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")



"""
#Print out sorted species file
finalnames = [Species_tag_dict[tag] for tag in Species_tag_dict]
final_cleaned_df = pd.DataFrame(columns=finalnames)
"""

"""
FUNCTIONS
"""

"""
Create VCF for snps dict
"""
"""
VCF_snps_dict = {}
with gzip.open(ind_snps, "rt") as in_snps:
    for line in in_snps:
        if not line.startswith("#"):
            fields = line.rstrip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            ident = chrom + "_" + pos
            ref = fields[3]
            alt = fields[4]
            genotype = fields[9].split(":")[0]
            VCF_snps_dict[ident] = [ref, alt, genotype]
"""
"""
Create VCF for indels dict


VCF_indels_dict = {}
with gzip.open(ind_snps, "rt") as in_snps:
    for line in in_snps:
        if not line.startswith("#"):
            fields = line.rstrip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            ident = chrom + "_" + pos
            ref = fields[3]
            alt = fields[4]
            genotype = fields[9].split(":")[0]
            VCF_indels_dict[ident] = [ref, alt, genotype]

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

"""
Function to obtain genotype from vcfs dicts
"""

def get_genotype(alleles_list):
    if alleles_list[2] == "0/0":
        gt = alleles_list[0]
    elif alleles_list[2] == "0/1" or alleles_list[2] == "1/0" :
        gt = alleles_list[0]
    elif alleles_list[2] == "1/1":
        gt = alleles_list[1]
    else:
        gt = alleles_list[0]
    return gt



"""
Parse fasta into dictionary
"""
records = SeqIO.to_dict(SeqIO.parse(open(reference_fasta), 'fasta'))


#MAIN SCRIPT
total_gene_seq = ""
ref_gene_seq = ""
with open(reference_gff, "r") as in_gff, \
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
