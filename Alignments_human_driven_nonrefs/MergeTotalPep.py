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
    merged_fasta_path = sys.argv[1]
    refs_species = sys.argv[2]
    out_file = sys.argv[3]
    duplicated_species = sys.argv[4]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

#DICT SPECIES TAGS REFERENCES
target_species = refs_species['Species'].to_list()
tag_species = refs_species['Tag'].to_list()
Species_tag_dict = {tag_species[i]:target_species[i] for i in range(0, len(tag_species))}

#dict DUPLICATED species
Duplicated_df = pd.read_csv(duplicated_species, sep='\t', low_memory=False)#panda creation
species_to_exclude = Duplicated_df['Species'].to_list()
tags = Duplicated_df['Ref_tag'].to_list()
Duplicated_tag_dict = {species_to_exclude[i]:tags[i] for i in range(0, len(species_to_exclude))}



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
species_idents = []
species_to_print = False
"""
let's remove duplicates species between references and resequenced data
"""
with gzip.open(merged_fasta_path, 'rt') as f, gzip.open(final_fasta_path, 'wt') as out_fh:
    for line in f:
        line = line.rstrip()
        if line.startswith(">") and any(line[1:].startswith(tag_id) for tag_id in tag_species):
            ref_species = select_current_value(Species_tag_dict, line[1:])
            print(">"+ref_species, file=out_fh)
            species_to_print = True
            ref_species_idents.append(ident)
        elif line.startswith(">") and not any(line[1:].startswith(spc) for spc in species_to_exclude):
            print(line, file=out_fh)
            species_to_print = True
        elif line.startswith(">") and any(line[1:].startswith(spc) for spc in species_to_exclude):
            corresponding_tag = select_current_value(Duplicated_tag_dict, line[1:])
            if any(ident.startswith(corresponding_tag) for ident in ref_species_idents):
                species_to_print = False
                continue
            else:
                print(line, file=out_fh)
                species_to_print = True
        elif not line.startswith(">") and species_to_print:
            print(line, file=out_fh)
