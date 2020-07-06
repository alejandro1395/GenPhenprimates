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

if len(sys.argv) == 9:
    alignment_path = sys.argv[1]
    pep_fasta_path =  sys.argv[2]
    pep_resequenced_path = sys.argv[3]
    tags_file = sys.argv[4]
    cds_fasta_path = sys.argv[5]
    cds_resequenced_path = sys.argv[6]
    out_path_cds = sys.argv[7]
    out_path_pep = sys.argv[8]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

#SPECIES NAMES ASSOCIATES TO TAGS
column1 = "Tag"
column2 = "SpeciesBroad"
Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
Ref_tags = pd.read_csv(cds_fasta_path+"summary_species.txt", sep='\t', low_memory=False)

#DICT SPECIES TAGS
Genomes_names = Species_tags['Genomes_names'].to_list()
Tree_names = Species_tags['Tree_names'].to_list()
Species_tag_dict = {Genomes_names[i]:Tree_names[i] for i in range(0, len(Genomes_names))}
Reference_species = Ref_tags['Species'].to_list()
Reference_idents = Ref_tags['Tag'].to_list()
Refs_tag_dict = {Reference_species[i]:Reference_idents[i] for i in range(0, len(Reference_species))}
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

"""
Function to check if reference in species fasta or take from resequenced data
"""

def check_reference_in_species_fasta(read_path, value):
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
    return target_sequence



#MAIN
#Create FASTA FILES
used_species = []
with open(alignment_path, "rt") as in_fh, gzip.open(out_path_cds, 'wt') as out_fh, \
gzip.open(out_path_pep, 'wt') as out_fh2:
    for line in in_fh:
        spc_line = line.rstrip()
        if (any(spc_line.startswith(str(genome)) for genome in Genomes_names)):
            curr_species = spc_line.split(" ")[0]
            """
            #REFERENCES GENOME ALGORITHM
            """
            if curr_species in Reference_species:
                pep_path = pep_fasta_path
                cds_path = cds_fasta_path + curr_species + "/" + curr_species + ".cds.fa.gz"
                if curr_species not in used_species:
                    tag_species = select_current_value(Refs_tag_dict, curr_species)
                    presence = check_reference_in_species_fasta(pep_fasta_path, tag_species)
                    if presence:
                        ident_pep, seq_pep = select_sequence_from_species_fasta(pep_path, tag_species)
                        print(">"+curr_species, file=out_fh2)
                        print(seq_pep, file=out_fh2)
                        id_spc = ident_pep.split(" ")[0][1:]
                        if id_spc.startswith("Propithecus_coquereli") or id_spc.startswith("Carlito_syrichta"):
                            id_spc = id_spc[:-2]
                        ident_cds, seq_cds = select_sequence_from_species_fasta(cds_path, id_spc)
                        print(">"+curr_species, file=out_fh)
                        print(seq_cds, file=out_fh)
            else:
                """
                #CASE OF RESEQUENCED SPECIES
                """
                pep_path = pep_resequenced_path
                cds_path = cds_resequenced_path
                if curr_species not in used_species:
                    ident_pep, seq_pep = select_sequence_from_species_fasta(pep_path, curr_species)
                    print(">"+curr_species, file=out_fh2)
                    print(seq_pep, file=out_fh2)
                    ident_cds, seq_cds = select_sequence_from_species_fasta(cds_path, curr_species)
                    print(">"+curr_species, file=out_fh)
                    print(seq_cds, file=out_fh)
            used_species.append(curr_species)
