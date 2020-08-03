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
import tarfile
import glob

"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 7:
    references_path = sys.argv[1]
    reference_names = sys.argv[2]
    tree_names = sys.argv[3]
    resequenced_path = sys.argv[4]
    out_path = sys.argv[5]
    gene_name = sys.argv[6]
else:
    sys.exit("The usage shoud be: human_path primates_path out_path")

"""
#VARIABLES
"""

#CREATE DICTS OF NAMES FOR REFERENCES
Species_tags = pd.read_csv(reference_names, sep='\t', low_memory=False)#panda creation
target_species = Species_tags['Species'].to_list()
tag_species = Species_tags['Tag'].to_list()
References_tag_dict = {tag_species[i]:target_species[i] for i in range(0, len(tag_species))}

#CREATE DICTS OF NAMES OF TREE
tree_tags = pd.read_csv(tree_names, sep='\t', low_memory=False)#panda creation
Genomes_names = tree_tags['Genomes_names'].to_list()
Tree_names = tree_tags['Tree_names'].to_list()
Species_tag_dict = {Tree_names[i]:Genomes_names[i] for i in range(0, len(Genomes_names))}

"""

#SPECIES NAMES ASSOCIATES TO TAGS
Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames

#DICT SPECIES IDENTITIES
Species_ident_dict = {}
Species_mean_dict = {}
target_species = Species_tags['Species'].to_list()
for i in range(0, len(target_species)):
    Species_ident_dict[target_species[i]] = {}
    Species_mean_dict[target_species[i]] = []
    for j in range(0, len(target_species)):
        Species_ident_dict[target_species[i]][target_species[j]] = []

Genes_mean_ident = []
Species_mean_dict["Gene"] = []

#DICT TAG IDENTITIES
tag_species = Species_tags['Tag'].to_list()
Species_tag_dict = {tag_species[i]:target_species[i] for i in range(0, len(tag_species))}

#all primates list
all_primates_dfs = {}


FUNCTIONS


def isNaN(string):
    return string != string


#Create gene families context with humans
human_families = pd.read_csv(human_path, sep='\t', low_memory=False)
human_genes = {}
print(human_families)

for index, row in human_families.iterrows():
    if row['hgnc_symbol'] not in human_genes and not isNaN(row['hgnc_symbol']):
        human_genes[row['hgnc_symbol']] = {}
        human_genes[row['hgnc_symbol']]["Homo_sapiens"] = []
        human_genes[row['hgnc_symbol']]["Homo_sapiens"].append(row['ensembl_peptide_id'])
    elif not isNaN(row['hgnc_symbol']):
        human_genes[row['hgnc_symbol']]["Homo_sapiens"].append(row['ensembl_peptide_id'])

#print(human_genes)
"""


#print(human_genes)

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


#tar files to avoid multiple gene problem
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


def select_current_value(current_dict, name):
    sel_value = [current_dict[str(tag)] for tag in current_dict if name.startswith(str(tag))][0]
    return sel_value



#MAIN
#Create FASTA FILES
#Now print the created files

ref_sequences = {}
with gzip.open(references_path, "rt") as in_fh:
    for line in in_fh:
        line = line.rstrip()
        if line.startswith(">"):
            ident = line
            ref_sequences[ident[1:]] = ""
        else:
            ref_sequences[ident[1:]] += line


#PRINT REFERENCE TOGETHER WITH THEIR INDIVIDUALS
for ref in ref_sequences:
    ref_genome = Species_tag_dict[ref]
    for k, v in References_tag_dict.items():
        if v == ref_genome:
            value = k
    presence = check_reference_in_species_fasta(resequenced_path, value)
    if presence:
        with open(out_path+ref+".within.cds.fa", "w") as out_fh:
            print(">"+ref, file=out_fh)
            print(ref_sequences[ref], file=out_fh)
            match_spc = False
            with open(resequenced_path, 'r') as f:
                for line in f:
                    line = line.rstrip()
                    if line.startswith(">"+value):
                        print(line, file=out_fh)
                        match_spc = True
                    elif line.startswith(">"):
                        match_spc = False
                    elif not line.startswith(">") and match_spc:
                        print(line, file=out_fh)
make_tarfile(out_path+gene_name+".tar.gz", out_path)
os.system('rm -rf '+out_path+'*.cds.fa')
