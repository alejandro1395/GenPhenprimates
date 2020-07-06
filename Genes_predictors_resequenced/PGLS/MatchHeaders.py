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

if len(sys.argv) == 6:
    input_path_cds = sys.argv[1]
    input_path_pep = sys.argv[2]
    tree_names_path = sys.argv[3]
    output_path_cds = sys.argv[4]
    output_path_pep = sys.argv[5]

else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

#SPECIES NAMES ASSOCIATES TO TAGS
column1 = "Tag"
column2 = "SpeciesBroad"
Tree_names = pd.read_csv(tree_names_path, sep='\t', low_memory=False)#panda creation

#DICT SPECIES TAGS
Genomes_names = Tree_names['Genomes_names'].to_list()
Tree_names = Tree_names['Tree_names'].to_list()
Tree_names_dict = {Genomes_names[i]:Tree_names[i] for i in range(0, len(Tree_names))}


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

def select_current_value(current_dict, name):
    sel_value = [current_dict[tag] for tag in current_dict if name.startswith(str(tag))][0]
    return sel_value

"""
Function to select fasta from species file and store it in a dict format
"""

def print_new_fasta_names_from_tree(read_path, out_path):
    gene_dict = {}
    skip = False
    with gzip.open(read_path, 'rt') as f, gzip.open(out_path, "wt") as out_f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                genome_species = line[1:]
                tree_species = select_current_value(Tree_names_dict, genome_species)
                if isNaN(tree_species):
                    skip = True
                    continue
                else:
                    ident = ">"+tree_species
                    print(ident, file=out_f)
                    skip = False
            elif not line.startswith(">") and not skip:
                print(line, file=out_f)


#MAIN
#Create FASTA FILES
print_new_fasta_names_from_tree(input_path_cds, output_path_cds)

print_new_fasta_names_from_tree(input_path_pep, output_path_pep)
