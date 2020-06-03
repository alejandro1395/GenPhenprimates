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


"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 7:
    alignment_path = sys.argv[1]
    pep_fasta_path = sys.argv[2]
    tree_file = sys.argv[3]
    cds_fasta_path = sys.argv[4]
    out_path_cds = sys.argv[5]
    out_path_pep = sys.argv[6]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

#SPECIES NAMES ASSOCIATES TO TAGS
column1 = "Tag"
column2 = "SpeciesBroad"
Tree_names = pd.read_csv(tree_file, sep='\t', low_memory=False)#panda creation


#DICT SPECIES TAGS
Genomes_names = Tree_names['Genomes_names'].to_list()
Tree_names = Tree_names['Tree_names'].to_list()
Tree_names_dict = {Tree_names[i]:Genomes_names[i] for i in range(0, len(Tree_names))}
"""
#Print out sorted species file
finalnames = [Species_tag_dict[tag] for tag in Species_tag_dict]
final_cleaned_df = pd.DataFrame(columns=finalnames)
"""

"""
FUNCTIONS
"""

"""
Function to check nan values
"""

def isNaN(string):
    return string != string

"""
Select current species from all species dataframes
"""

def select_current_value(current_dict, name, col):
    sel_value = [current_dict[tag] for tag in current_dict if name.startswith(tag)][0]
    return sel_value

"""
Function to select fasta from species file and store it in a dict format
"""

def select_sequence_from_species_fasta(read_path, value):
    spc_dict = {}
    with gzip.open(read_path, 'rt') as f:
        target_sequence = False
        for line in f:
            line = line.rstrip()
            if line.startswith(">"+value):
                ident = line
                spc_dict[ident] = ""
                target_sequence = True
            elif line.startswith(">") and target_sequence:
                break
            elif not line.startswith(">") and target_sequence:
                spc_dict[ident] += line
    return ident, spc_dict[ident]


#MAIN
#Create FASTA FILES
used_species = []
with open(alignment_path, "rt") as in_fh, gzip.open(out_path_cds, 'wt') as out_fh, \
gzip.open(out_path_pep, 'wt') as out_fh2:
    for line in in_fh:
        spc_line = line.rstrip()
        if (any(spc_line.startswith(tree_name) for tree_name in Tree_names_dict)):
            new_fields = spc_line.split(" ")
            curr_species = new_fields[0]
            pep_path = pep_fasta_path
            cds_path = cds_fasta_path
            if curr_species not in used_species:
                ident_pep, seq_pep = select_sequence_from_species_fasta(pep_path, new_fields[0])
                print(ident_pep, file=out_fh2)
                print(seq_pep, file=out_fh2)
                id_spc = ident_pep.split(" ")[0][1:]
                if id_spc.startswith("Propithecus_coquereli") or id_spc.startswith("Carlito_syrichta"):
                    id_spc = id_spc[:-2]
                print(id_spc)
                ident_cds, seq_cds = select_sequence_from_species_fasta(cds_path, id_spc)
                print(ident_cds)
                print(ident_cds, file=out_fh)
                print(seq_cds, file=out_fh)
            used_species.append(curr_species)
