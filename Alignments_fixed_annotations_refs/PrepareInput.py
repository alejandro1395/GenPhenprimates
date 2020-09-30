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

if len(sys.argv) == 5:
    clusters_path = sys.argv[1]
    tags_file = sys.argv[2]
    pep_fasta_path = sys.argv[3]
    out_path = sys.argv[4]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

#VARIABLES
raw_df = pd.read_csv(clusters_path, sep='\t', low_memory=False)#panda creation
column1 = "Tag"
column2 = "SpeciesBroad"


#SPECIES NAMES ASSOCIATES TO TAGS
Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames

#DICT SPECIES TAGS
target_species = Species_tags['Species'].to_list()
tag_species = Species_tags['Tag'].to_list()
Species_tag_dict = {tag_species[i]:target_species[i] for i in range(0, len(tag_species))}
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

def select_sequence_from_species_fasta(value, spc_dict):
    curr_species = select_current_value(Species_tag_dict, value, column1)
    read_path = pep_fasta_path + curr_species + "/" + curr_species + ".pep.fa.gz"
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
fastas_count = 0
for index, row in raw_df.iterrows():
    fastas_count += 1
    with gzip.open(out_path, 'wt') as f:
        #FUNCTION TO CREATE ALIGNMENT FILE
        for value in row.values.tolist():
            species_fasta = {}
            if isNaN(value):
                continue
            else:
                ID, sequence = select_sequence_from_species_fasta(value, species_fasta)
                print(ID, file=f)
                print(sequence, file=f)
