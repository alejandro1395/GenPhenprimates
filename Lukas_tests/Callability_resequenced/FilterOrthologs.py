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
    all_genes = sys.argv[1]
    ortholog_list = sys.argv[2]
    species_name_list = sys.argv[3]
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
Read original files
"""

gene_overlap_info = pd.read_csv(all_genes, sep='\t', low_memory=False,
names=["Individual", "Gene", "Length", "Prop"])

ortholog_list_df = pd.read_csv(ortholog_list, sep='\t', low_memory=False,
names=["Gene"])
all_orthologs = ortholog_list_df['Gene'].to_list()


def select_current_value(current_dict, name):
    sel_value = [current_dict[tag] for tag in current_dict if name.startswith(tag)][0]
    return sel_value

"""
SPECIES NAME gene_overlap_info
"""

Species_tags = pd.read_csv(species_name_list, sep='\t', low_memory=False)
#DICT SPECIES TAGS
individuals_ids = Species_tags['PDGP_ID'].to_list()
Species_tags['Species_name'] = Species_tags[['GENUS', 'SPECIES']].agg('_'.join, axis=1)
species_names = Species_tags['Species_name'].to_list()
Species_tag_dict = {individuals_ids[i]:species_names[i] for i in range(0, len(individuals_ids))}


#Now let's check if the ortho is in the list,
#then print it

ortho_overlap_final_df = pd.DataFrame(columns=("Individual", "Gene", "Length", "Prop", "Species"))

for i, row in gene_overlap_info.iterrows():
    gene_id = str(row['Gene'])
    if gene_id in all_orthologs:
        species = select_current_value(Species_tag_dict, row['Individual'])
        values_to_add = {'Individual':  row['Individual'], 'Gene': gene_id,
        'Length': row['Length'], 'Prop': row['Prop'], 'Species': species}
        row_to_add = pd.Series(values_to_add)
        ortho_overlap_final_df = ortho_overlap_final_df.append(row_to_add, ignore_index = True)

ortho_overlap_final_df.to_csv(out_file, sep="\t")
