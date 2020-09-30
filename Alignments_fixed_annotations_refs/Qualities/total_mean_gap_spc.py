#!/usr/bin/python3

import pandas as pd
import numpy as np
from numpy import array
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
Genes_cumulative_length = []
Genes_cumulative_gap = []

#DICT TAG IDENTITIES
tag_species = Species_tags['Tag'].to_list()
target_species = Species_tags['Species'].to_list()
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
alignments_count = 0
species = [e for e in dirs if e not in ('qu', 'out')]
for ref in species:
    for file in os.listdir(trimals_path+"/"+ref):
        if file.endswith("SpeciesGaps"):
            with open(trimals_path+"/"+ref+"/"+file, 'r') as in_fh:
                species_count = 0
                start_info = False
                for line in in_fh:
                    line = line.rstrip()
                    fields = line.split("\t")
                    if start_info:
                        new_fields = line.split('\t')
                        Genes_cumulative_length.append(float(new_fields[6]))
                        Genes_cumulative_gap.append(float(new_fields[10]))
                    if line.startswith("+"):
                        start_info = True

"""
#PRINT WHOLE VALUES FOR DISTRIBUTIONS
distributions = pd.DataFrame(Species_mean_dict)
distributions.to_csv(out_path+"gene_distr_mean_ident_spc", sep = "\t", index=False)
"""

#PRINT WHOLE VALUES FOR GENE DISTRIBUTIONS#create new df
distributions_gene = pd.DataFrame({"Genes_cumulative_length":Genes_cumulative_length,
"Genes_cumulative_gap":Genes_cumulative_gap})
distributions_gene.to_csv(out_path+"gene_distr_mean_gap", sep = "\t", index=False)

"""
#CALCULATE MEANS FOR ALL PAIR-SPECIES IDENTITIES
for spc1 in Species_ident_dict:
    for spc2 in Species_ident_dict[spc1]:
        Species_ident_dict[spc1][spc2] = mean(Species_ident_dict[spc1][spc2])

#BUILD MATRIX
matrix_species = pd.DataFrame(Species_ident_dict)
matrix_species.to_csv(out_path+"total_mean_ident_spc.tsv", sep = "\t", index=False)
"""
