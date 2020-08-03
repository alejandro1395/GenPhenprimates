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
    gene_path = sys.argv[1]
    gene_name = sys.argv[2]
    out_path_pep = sys.argv[3]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

#SPECIES NAMES ASSOCIATES TO TAGS
column1 = "Tag"
column2 = "SpeciesBroad"
#Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation

#DICT SPECIES TAGS
#DICT SPECIES TAGS
"""
individuals_ids = Species_tags['PGDP ID'].to_list()
Species_tags['Species_name'] = Species_tags[['Genus', 'Species']].agg('_'.join, axis=1)
species_names = Species_tags['Species_name'].to_list()
Species_tag_dict = {individuals_ids[i]:species_names[i] for i in range(0, len(individuals_ids))}
"""
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
    sel_value = [current_dict[tag] for tag in current_dict if name.startswith(tag)][0]
    return sel_value

"""
Function to select fasta from species and count intermediate stop codons
"""

def count_stops_in_species_fasta(read_path):
    spc_dict = {}
    with open(read_path, 'r') as f:
        target_sequence = False
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                ident = line[1:]
                spc_dict[ident] = 0
            elif not line.startswith(">"):
                for element in line:
                    if element == "*":
                        spc_dict[ident] += 1
    return spc_dict

"""
Function to count percentage of protein without stop codons
"""

def count_prot_uninterrupted_in_species_fasta(read_path):
    spc_dict_seq = {}
    spc_dict_perc = {}
    ident = ""
    with open(read_path, 'r') as f:
        target_sequence = False
        for line in f:
            line = line.rstrip()
            if line.startswith(">") and ident:
                total = len(spc_dict_seq[ident])
                if "*" in spc_dict_seq[ident]:
                    first_stop_codon = spc_dict_seq[ident].find("*")
                    good_perc = (int(first_stop_codon)+1)/total
                    spc_dict_perc[ident] = good_perc
                else:
                    spc_dict_perc[ident] = 1
                ident = line[1:]
                spc_dict_seq[ident] = ""
            elif line.startswith(">"):
                ident = line[1:]
                spc_dict_seq[ident] = ""
            elif not line.startswith(">"):
                spc_dict_seq[ident] += line
        total = len(spc_dict_seq[ident])
        if "*" in spc_dict_seq[ident]:
            first_stop_codon = spc_dict_seq[ident].find("*")
            good_perc = (int(first_stop_codon)+1)/total
            spc_dict_perc[ident] = good_perc
        else:
            spc_dict_perc[ident] = 1
    return spc_dict_perc




#MAIN
#Create FASTA FILES
species_df = pd.DataFrame(columns=['Gene', 'Species', 'Individual', 'Count', 'Good_perc'])
individuals_stops = count_stops_in_species_fasta(gene_path)
individuals_perc = count_prot_uninterrupted_in_species_fasta(gene_path)
for ind in individuals_stops:
    species_df = species_df.append({'Gene': gene_name, 'Species': ind.split("_")[0], 'Individual': ind,
            'Count': individuals_stops[ind], 'Good_perc': individuals_perc[ind]}, ignore_index = True)

output_file = out_path_pep
species_df.to_csv(output_file, sep='\t')
