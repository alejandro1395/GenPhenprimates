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
import fileinput


"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 4:
    tags_file = sys.argv[1]
    ref_species = sys.argv[2]
    out_file = sys.argv[3]
else:
    sys.exit("The usage shoud be: ./FinderSBH.py in_file tag_file output_file")

#VARIABLES
#!/usr/bin/env python3

#VARIABLES
query_species = ""
SBH_dict = {}
species_counter = {}
cluster_dict = {}

"""
Read TAGs file foreach primate species used in the analysis
and store it in a Pandas DataFrame
"""

Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames
SBH_df = pd.DataFrame(columns=finalnames)


"""
FUNCTIONS
"""

"""
Store in a dictionary all target species so as to keep a counter
for unique species as best hits
"""

def store_target_species_count_in_dict(Species_df, reference):
    target_species = Species_df['Species'].to_list()
    print(reference)
    target_species.remove(reference)
    species_counter = {name:0 for name in target_species}
    return species_counter


"""
Here, we create a function that parses the clusters of sequences to retrieve all
IDs from the cluster appart from the representative one
"""

def parse_cluster_file(clusters_file, clusters_dict):
    with open(clusters_file, "rt") as in_fh:
        id = ""
        clustered = []
        for line in in_fh:
            line = line.rstrip()
            if line.startswith(">") and not id:
                continue
            elif line.startswith(">"):
                clusters_dict[representative] = clustered
                id = ""
                clustered = []
            else:
                id = line.split(">")[1].split(".")[0]
                if line[-1] == "*":
                    representative = id
                else:
                    clustered.append(id)
        clusters_dict[representative] = clustered
    return clusters_dict


"""
Another function here helps us to identify extra species by looking whether the
ID of the species match any of the clusters, and retrieves all the other
target IDs associated
"""


def check_all_species_in_cluster(ident, query, clusters_dict, best_hit_dict,
sp_counter):
    sp_counter, best_hit_dict = include_only_best_hit_foreach_species_target(ident,
    best_hit_dict, query, sp_counter)
    if ident in clusters_dict:
        for value in clusters_dict[ident]:
            sp_counter, best_hit_dict = include_only_best_hit_foreach_species_target(value,
            best_hit_dict, query, sp_counter)
    return sp_counter, best_hit_dict


"""
Here, the function includes only a species once, as the best hit for our reference
species. Then, with that in mind, filters for species appearing more than once
as hits for a specific query entry
"""

#FUNCTION TO INCLUDE ONLY A SPECIES ONCE (BEST HIT)
def include_only_best_hit_foreach_species_target(elem, best_hit_dict, query,
sp_counter):
    for item in Species_tags['Tag']:
        if elem.startswith(item):
            current_species = Species_tags.loc[Species_tags['Tag'] == item].Species.item()
            if current_species in sp_counter:
                if sp_counter[current_species] == 0:
                    sp_counter[current_species] += 1
                    best_hit_dict[query].append(elem)
    return sp_counter, best_hit_dict


"""
Last function to print the ouput in Pandas format for the Query and
Target columns in our dataframe
"""

#FUNCTION TO INCLUDE ONLY A SPECIES ONCE (BEST HIT)
def append_out_BBHs_pandas_format(sbh_dict, sbh_df, query):
    query_row = [query] + sbh_dict[query] + \
    list(np.repeat(np.nan, len(sbh_df.columns)-len(sbh_dict[query])-1))
    sbh_df = sbh_df.append(pd.Series(query_row, index=sbh_df.columns),
    ignore_index=True)
    return sbh_df



"""
MAIN
"""

species_counter = store_target_species_count_in_dict(Species_tags, ref_species)

#cluster_dict = parse_cluster_file(clust_file, cluster_dict)

count = 0
in_fh = iter(sys.stdin)
for line in in_fh:
    line = line.rstrip()
    if line.startswith("Query=") and query_species == "":
        query_fields = line[7:].split(" ")
        query_species = "_".join(query_fields[0:1])
        SBH_dict[query_species] = []
    elif line.startswith("Query="):
        SBH_df = append_out_BBHs_pandas_format(SBH_dict, SBH_df,
        query_species)
        query_fields = line[7:].split(" ")
        query_species = "_".join(query_fields[0:1])
        SBH_dict[query_species] = []
        species_counter = {k:0 for k in species_counter}
    elif line.startswith("Sequences producing significant alignments"):
        next(in_fh)
        line_new = next(in_fh).rstrip()
        while (any(letter.isalnum() for letter in line_new)):
            ID_fields = line_new.split(" ")
            ID = "_".join(ID_fields[2:3])
            species_counter, SBH_dict = include_only_best_hit_foreach_species_target(ID, SBH_dict,
            query_species, species_counter)
            line_new = next(in_fh).rstrip()

#REPEAT THIS AFTER LOOP`FOR LAST HOMOLOG ENTRY
SBH_df = append_out_BBHs_pandas_format(SBH_dict, SBH_df,
query_species)
#Print_SBHs_in_Pandas_format(SBH_dict, SBH_df, out_file)
SBH_df.to_csv(out_file, sep = "\t", index=False)
