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

if len(sys.argv) == 5:
    human_genes = sys.argv[1]
    reference_gff = sys.argv[2]
    coord_high_qual = sys.argv[3]
    out_file = sys.argv[4]
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
Read human genes
"""

with open(human_genes, "r") as in_fh:
    all_human_genes = in_fh.read().splitlines()

"""
Read high quality intervals
"""

qual_intervals = []
with open(coord_high_qual, "r") as in_fh2:
    for line in in_fh2:
        fields = line.rstrip().split("\t")
        qual_intervals.append([fields[0], fields[1], fields[2]])


"""
Function to detect overlapping intervals
"""

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


#MAIN SCRIPT


with gzip.open(reference_gff, "rt") as in_gff:
    last_line = in_gff.readlines()[-1]

genes_dict = {}
last_Parent_info = None
with gzip.open(reference_gff, "rt") as in_gff:
    for line in in_gff:
        if line == last_line:
            break
        fields = line.rstrip().split("\t")
        line = line.rstrip()
        chrom = fields[0]
        Parent_info = fields[8].split("=")[1].split(";")[0]
        seq_type = fields[2]
        if seq_type == "CDS" and Parent_info in all_human_genes:
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            if Parent_info not in genes_dict and last_Parent_info is not None:
                genes_dict[last_Parent_info] = genes_dict[last_Parent_info]/total_length
                total_length = 0
                genes_dict[Parent_info] = 0
            elif Parent_info not in genes_dict:
                total_length = 0
                genes_dict[Parent_info] = 0
            curr_interval = [chrom, start, end]
            seq_length = end - start
            total_length += seq_length
            interval_found = False
            for inter in qual_intervals:
                if curr_interval[0] == "chr"+inter[0]:
                    overlap = getOverlap([int(curr_interval[1]), int(curr_interval[2])],
                    [int(inter[1]), int(inter[2])])
                    if overlap != 0:
                        if (int(curr_interval[1]) >= int(inter[1])) and (int(curr_interval[2]) <= int(inter[2])):
                            subset = int(curr_interval[2]) - int(curr_interval[1])
                        elif (int(curr_interval[1]) <= int(inter[1])) and (int(curr_interval[2]) <= int(inter[2])):
                            subset = int(curr_interval[2]) - int(inter[1])
                        elif (int(curr_interval[1]) <= int(inter[1])) and (int(curr_interval[2]) >= int(inter[2])):
                            subset = int(inter[2]) - int(inter[1])
                        else:
                            subset = int(inter[2]) - int(curr_interval[1])
                        genes_dict[Parent_info] += subset
            last_Parent_info = Parent_info


Genes = pd.DataFrame(columns=("Gene", "Perc"))
for gene in genes_dict:
    values_to_add = {'Gene': gene, 'Perc': round(genes_dict[gene], 2)}
    row_to_add = pd.Series(values_to_add)
    Genes = Genes.append(row_to_add, ignore_index = True)
#Spc_frame = pd.DataFrame.from_dict(SpeciesRefs_frameshift_dict, orient='index')
Genes.to_csv(out_file, sep="\t")
