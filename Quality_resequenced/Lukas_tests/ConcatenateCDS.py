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
import argparse
import glob



"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 5:
    cds_fasta = sys.argv[1]
    reference_gff = sys.argv[2]
    gene_name = sys.argv[3]
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
Parse gff into dictionary to detect strand
"""

CDS_coord_strand = {}
with open(reference_gff, "r") as in_gff:
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
            print(chrom+":"+str(start)+"-"+str(end))
            CDS_coord_strand[chrom+":"+str(start)+"-"+str(end)] = strand


"""
Sort numerically the chunks by number
"""

def SortChunks(value_list):
    sorted_new_list = []
    for i in range(1, len(value_list)+1):
        for val in value_list:
            pattern = "chunk" + str(i) + "."
            if pattern in val:
                sorted_new_list.append(val)
    return sorted_new_list

"""
Parse all fasta files to concatenate
"""
list_files = glob.glob(cds_fasta+'/'+'*.fa')
sorted_list = SortChunks(list_files)

sequence = ""
for file in sorted_list:
    with open(file, "r") as in_fh:
        for line in in_fh:
            line = line.rstrip()
            if line.startswith(">"):
                curr_strand = CDS_coord_strand[line[1:]]
            else:
                sequence += line

with open(out_file+"/"+gene_name+".fa", "w") as out_fh:
    if curr_strand == "-":
        final_seq = reverse_complement(sequence)
    else:
        final_seq = sequence
    print(">"+gene_name, file=out_fh)
    print(final_seq, file=out_fh)
    exit()


#rc_con_seq = Seq(str(con_seq)).reverse_complement()
#print rc_con_seq, len(rc_con_seq)
