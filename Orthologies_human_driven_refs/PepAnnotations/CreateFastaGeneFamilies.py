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

if len(sys.argv) == 5:
    human_path = sys.argv[1]
    tags_file = sys.argv[2]
    primates_path = sys.argv[3]
    out_path = sys.argv[4]
else:
    sys.exit("The usage shoud be: human_path primates_path out_path")

#VARIABLES
column1 = "Tag"
column2 = "SpeciesBroad"


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

"""
FUNCTIONS
"""

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

#Add primates info to these gene families
dirs = os.listdir(primates_path)
species = [e for e in dirs if e not in ('qu', 'out', "summary_species.txt")]
for spc in species:
    for file in os.listdir(primates_path+spc):
        if file.endswith("ref_genes.tsv"):
            path_to_pandas = primates_path + spc + "/" + spc + ".ref_genes.tsv"
            df_species = pd.read_csv(path_to_pandas, index_col=None,
            sep ="\t", low_memory=False, header=None, names=["GeneID", "Gene_name"])
            all_primates_dfs[spc] = df_species

#print(all_primates_dfs)

for primate in all_primates_dfs:
    if primate != "Homo_sapiens":
        for index, row in all_primates_dfs[primate].iterrows():
            if row['Gene_name'] in human_genes and primate not in human_genes[row['Gene_name']]:
                human_genes[row['Gene_name']][primate] = []
                human_genes[row['Gene_name']][primate].append(row['GeneID'])
            elif row['Gene_name'] in human_genes and primate in human_genes[row['Gene_name']]:
                human_genes[row['Gene_name']][primate].append(row['GeneID'])

#print(human_genes)

#Now recover the fastas from this dictionary
#Function to select fasta from species file and store it in a dict format


def select_sequence_from_species_fasta(read_path, ident, spc_dict, spc):
    if spc in ["Propithecus_coquereli", "Carlito_syrichta"]:
        ident = ident+"_1"
    identifier = ""
    with gzip.open(read_path, 'rt') as f:
        target_sequence = False
        for line in f:
            line = line.rstrip()
            fields = line.split(" ")
            if line.startswith(">") and fields[0].endswith(ident):
                identifier = line
                spc_dict[identifier] = ""
                target_sequence = True
            elif line.startswith(">") and target_sequence:
                break
            elif not line.startswith(">") and target_sequence:
                spc_dict[identifier] += line
        if not identifier:
            identifier = "null"
            spc_dict[identifier] = "null"
    return identifier, spc_dict[identifier]

#tar files to avoid multiple gene problem
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

#MAIN
#Create FASTA FILES
#Now print the created files

for gene in human_genes:
    if not os.path.exists(out_path+gene):
        os.mkdir(out_path+gene)
    for primate in human_genes[gene]:
        if not os.path.exists(out_path+gene+"/"+primate):
            os.mkdir(out_path+gene+"/"+primate)
        for file in os.listdir(primates_path+primate):
            if file.endswith("pep.fa.gz"):
                with gzip.open(out_path+gene+"/"+primate+"/"+ primate+".pep.fa.gz", 'wt') as out_fh:
                    spc_gene_dict = {}
                    for value in human_genes[gene][primate]:
                        gene_id, sequence = select_sequence_from_species_fasta(primates_path+primate+"/"+file, value,
                        spc_gene_dict, primate)
                        if gene_id != "null" and sequence != "null":
                            print(gene_id, file=out_fh)
                            print(sequence, file=out_fh)
    make_tarfile(out_path+gene+"/"+gene+".tar.gz", out_path+gene+"/")
    os.system('rm -rf '+out_path+gene+'/*/')
