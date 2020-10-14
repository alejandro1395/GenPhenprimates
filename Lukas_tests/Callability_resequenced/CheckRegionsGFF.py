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
    ref_genes = sys.argv[1]
    individual = sys.argv[2]
    individual_overlap_gff = sys.argv[3]
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
Read lengths for original CDS files
"""
total_len_genes = {}
dirs = os.listdir(ref_genes)
species = [e for e in dirs if e not in ('qu', 'out', "summary_species.txt")]
for spc in species:
    for file in os.listdir(ref_genes+spc):
        if file.endswith("cds.fa.gz"):
            with gzip.open(ref_genes+spc+"/"+file, 'rt') as in_fh:
                for line in in_fh:
                    line = line.rstrip()
                    if line.startswith(">"):
                        ident = line[1:].split(" ")[0]
                        total_len_genes[ident] = 0
                    else:
                        total_len_genes[ident] += len(line)

"""
Read overlaps file
"""

overlap_dict = {}
with open(individual_overlap_gff, "r") as in_fh2:
    for line in in_fh2:
        fields = line.rstrip().split("\t")
        feature = fields[5]
        if feature == "CDS":
            ident = fields[11].split(";")[0].split("=")[1]
            overlap = int(fields[12])
            if ident in overlap_dict:
                overlap_dict[ident] += overlap
            else:
                overlap_dict[ident] = overlap

"""
print output
"""

with open(out_file, "w") as out_fh:
    for key in overlap_dict:
        prop = overlap_dict[key]/total_len_genes[key]
        print("{}\t{}\t{}\t{}".format(individual, key, overlap_dict[key], prop), file=out_fh)
