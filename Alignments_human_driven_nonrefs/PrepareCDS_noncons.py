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
    cds_fasta_path = sys.argv[1]
    gene_name = sys.argv[2]
    current_species = sys.argv[3]
    consensus = sys.argv[4]
    out_file = sys.argv[5]
    to_exclude = sys.argv[6]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

"""
Read file for species to exclude
"""

Species_to_exclude = pd.read_csv(to_exclude, sep='\t', low_memory=False)#panda creation
individuals_to_exclude = Species_to_exclude['PGDP ID'].to_list()





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


#Select current species from all species dataframes


def select_current_value(current_dict, name):
    sel_value = [current_dict[tag] for tag in current_dict if name.startswith(tag)][0]
    return sel_value


#Function to select fasta from species file and store it in a dict format

cons_seq = ""
sequences = {}
good_ind = False
with gzip.open(cds_fasta_path, 'rt') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">") and len(sequences) == 1:
            break
        elif line.startswith(">") and not any(line[1:].startswith(ind) for ind in individuals_to_exclude):
            ident = line
            sequences[ident] = ""
            good_ind = True
        elif not line.startswith(">") and good_ind:
            sequences[ident] += line
    with gzip.open(out_file+gene_name+".resequenced.cds.fa.gz", 'at') as out_f:
        for key in sequences:
            print(">"+current_species, file=out_f)
            print(sequences[key], file=out_f)
