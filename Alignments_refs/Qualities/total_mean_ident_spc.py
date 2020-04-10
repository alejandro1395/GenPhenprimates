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

if len(sys.argv) == 4:
    trimals_path = sys.argv[1]
    tags_file = sys.argv[2]
    out_path = sys.argv[3]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

#VARIABLES
column1 = "Tag"
column2 = "SpeciesBroad"


#SPECIES NAMES ASSOCIATES TO TAGS
Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames

#DICT SPECIES IDENTITIES
Species_ident_dict = {}
target_species = Species_tags['Species'].to_list()
for i in range(0, len(target_species)):
    Species_ident_dict[target_species[i]] = {}
    for j in range(0, len(target_species)):
        Species_ident_dict[target_species[i]][target_species[j]] = []
#DICT TAG IDENTITIES
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


"Function to iterate over all trimal alignments files"
def concatenate_BBH_all_species_into_dict(path, all_species_dfs):
    dirs = os.listdir(trimal)
    species = [e for e in dirs if e not in ('qu', 'out')]
    for folder in species:
        path_to_pandas = path + folder + "/" + folder + ".pep.BBH.gz"
        df_species = pd.read_csv(path_to_pandas, index_col=None,
        sep ="\t", low_memory=False)
        all_species_dfs.append(df_species)
    return all_species_dfs



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
#READ TRIMAL FILES
dirs = os.listdir(trimals_path)
species = [e for e in dirs if e not in ('qu', 'out')]
for ref in species:
    for file in os.listdir(trimals_path+"/"+ref):
        if file.endswith("SpeciesIdentities"):
            with open(trimals_path+"/"+ref+"/"+file, 'r') as in_fh:
                species_count = 0
                for line in in_fh:
                    line = line.rstrip()
                    fields = line.split("\t")
                    if line.startswith("#Percentage of identity matrix:"):
                        line_new = next(in_fh).rstrip()
                        used_species = []
                        while (any(line_new.startswith(tag) for tag in tag_species)):
                            new_fields = line_new.split("\t")
                            curr_species = select_current_value(Species_tag_dict, new_fields[0], column1)
                            used_species.append(curr_species)
                            species_count+=1
                            if len(used_species) == 1:
                                pass
                            else:
                                #Here we store the identity information for eahc pair of species twice
                                for count in range(0, len(used_species)):
                                    Species_ident_dict[curr_species][used_species[count]].append(float(new_fields[count+1].rstrip()))
                                    Species_ident_dict[used_species[count]][curr_species].append(float3444(new_fields[count+1].rstrip()))
                            line_new = next(in_fh).rstrip()
                    elif line.startswith("#Percentage of identity with most similar sequence:"):
                        print(Species_ident_dict)
                        exit()
