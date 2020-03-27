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
    clusters_path = sys.argv[1]
    tags_file = sys.argv[2]
    gff_path = sys.argv[3]
    alias_path = sys.argv[4]
    out_file = sys.argv[5]
else:
    sys.exit("The usage shoud be: ./OrthoCLust.py clusters_path tags_file out_file")


#VARIABLES
clusters_df = pd.read_csv(clusters_path, sep='\t', low_memory=False)#panda creation
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

#Print out sorted species file
finalnames = ["chr", "start", "end", "ClusterID", "score", "strand"]
final_cleaned_df = pd.DataFrame(columns=finalnames)

####LUEGO PENSARE SI LA INCLUYO
#validating species_(with more time add new options to validate)
val_species=["Chlorocebus_sabaeus", "Pan_troglodytes",
"Otolemur_garnettii", "Pan_paniscus", "Papio_anubis",
"Pongo_abelii", "Homo_sapiens"]
for species in val_species:
    final_cleaned_df.to_csv(out_file + species + ".Clusters_loc.pep.gz", header=False, index=False, sep="\t")


names_gff = ["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
names_alias = ["alias", "chrom", "source"]

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
#read pandas gff in loop for Validating species
"""

def import_files_from_validating_species(path, validating_species, end, names):
    validate_dict = {}
    dirs = os.listdir(path)
    for dir in dirs:
        if dir in validating_species:
            df_species = pd.read_csv(path + "/" + dir + "/" + dir + end,
            sep ="\t", names=names, low_memory=False)
            validate_dict[dir] = df_species
    return validate_dict


"""
Select coordinates info from dataframe of gff species
"""

def select_gff_info(dict1, dict2, species, value):
    df_sel1 = dict1[species]
    subset = df_sel1[df_sel1["attribute"].str.contains("ID="+value)]
    if species == "Homo_sapiens":
        chr = str(subset["chrom"].values[0])
    else:
        alias_chr = str(subset["chrom"].values[0])
        chr = select_value_from_alias(dict2, species, alias_chr)
    start = int(subset["start"].values)
    end = int(subset["end"].values)
    score = str(subset["score"].values[0])
    strand = str(subset["strand"].values[0])
    return chr, start, end, score, strand

"""
Sselect_value_from_alias in validating species
"""

def select_value_from_alias(dict2, species, value):
    df_sel2 = dict2[species]
    subset2 = df_sel2[df_sel2["alias"].str.contains(str(value))]
    print(subset2)
    chr = str(subset2["chrom"].values[0])
    return chr


"""
Function to parse all pairs of values and print them in the
Query vs multiple targets format
"""

def append_out_row_pandas_format(trans_row, output_file):
    trans_df = pd.DataFrame.from_dict(trans_row)
    print(trans_df)
    trans_df.to_csv(output_file, mode='a',header=False, index=False, sep="\t")

    #Drop nan values for 60-1 species, with clusters missing targets
    #(directed previously to other clusters)


#MAIN

gffs = import_files_from_validating_species(gff_path, val_species, ".gff.gz", names_gff)

chroms = import_files_from_validating_species(alias_path, val_species, ".chromAlias.txt.gz",
names_alias)

#Reorder based on species names as columns
clusters_count = 0
clusters_row = {}
for index, row in clusters_df.iterrows():
    clusters_count +=1
    for value in row.values.tolist():
        if not isNaN(value):
            curr_species = select_current_value(Species_tag_dict, value, column1)
            if curr_species in val_species:
                print(curr_species)
                species_row = {}
                Chr, Start, End, Score, Strand = select_gff_info(gffs, chroms,
                curr_species, value)
                species_row["chr"] = [Chr]
                species_row["start"] = [Start]
                species_row["end"] = [End]
                species_row["ClusterID"] = ["Cluster" + str(clusters_count)]
                species_row["score"] = [Score]
                species_row["strand"] = [Strand]
                out = out_file + curr_species + ".Clusters_loc.pep.gz"
                append_out_row_pandas_format(species_row, out)
