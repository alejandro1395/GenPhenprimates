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
    individuals_path = sys.argv[1]
    alignment_path = sys.argv[2]
    reference_names = sys.argv[3]
    tree_names = sys.argv[4]
    out_file = sys.argv[5]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")


#CREATE DICTS OF NAMES FOR REFERENCES
Species_tags = pd.read_csv(reference_names, sep='\t', low_memory=False)#panda creation
target_species = Species_tags['Species'].to_list()
tag_species = Species_tags['Tag'].to_list()
References_tag_dict = {tag_species[i]:target_species[i] for i in range(0, len(tag_species))}

#CREATE DICTS OF NAMES FOR TREE NAMES
Tree_tags = pd.read_csv(tree_names, sep='\t', low_memory=False)#panda creation
Tree_species = Tree_tags['Tree_names'].to_list()
Genome_species = Tree_tags['Genomes_names'].to_list()
Tree_tag_dict = {Tree_species[i]:Genome_species[i] for i in range(0, len(Genome_species))}



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
    sel_value = [current_dict[str(tag)] for tag in current_dict if name.startswith(str(tag))][0]
    return sel_value


#Function to select fasta from species file and store it in a dict format

"""
Function to select fasta from species file and store it in a dict format
"""

def select_sequence_from_species_fasta(read_path, value):
    spc_dict = {}
    with open(read_path, 'r') as f:
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

"""
Function to check if reference in species fasta or take from resequenced data
"""

def check_reference_in_species_fasta(read_path, value):
    spc_dict = {}
    with open(read_path, 'r') as f:
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
    return target_sequence



#MAIN
#Create FASTA FILES
used_genomes = []
with open(alignment_path, "r") as in_fh, open(out_file, 'a') as out_fh:
    for line in in_fh:
        spc_line = line.rstrip()
        if (any(spc_line.startswith(str(spc)) for spc in tag_species)):
            new_fields = spc_line.split(" ")
            curr_species_genome = select_current_value(References_tag_dict, new_fields[0])
            if curr_species_genome not in used_genomes:
                for filename in os.listdir(individuals_path):
                    ind = filename.split(".")[0]
                    presence = check_reference_in_species_fasta(individuals_path + "/" + filename, new_fields[0])
                    if presence:
                        ident_pep, seq_cds = select_sequence_from_species_fasta(individuals_path + "/" + filename, new_fields[0])
                        id_spc = ident_pep[1:]
                        if id_spc.startswith("Propithecus_coquereli") or id_spc.startswith("Carlito_syrichta"):
                            id_spc = id_spc[:-2]
                        print(ident_pep+"_"+ind, file=out_fh)
                        print(seq_cds, file=out_fh)
                used_genomes.append(curr_species_genome)
