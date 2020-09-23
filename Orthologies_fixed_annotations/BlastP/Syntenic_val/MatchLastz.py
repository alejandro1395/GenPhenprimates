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
    in_path = sys.argv[1]
    tags_file = sys.argv[2]
    out_file = sys.argv[3]
else:
    sys.exit("The usage shoud be: in_path tags_file non_lifted_path out_file")



#VARIABLES

clusters_dict_lift = {}
clusters_dict_nonlift = {}
Lifted_dfs = {}
Non_lifted_dfs = {}
val_species = val_species=["Chlorocebus_sabaeus", "Pan_troglodytes",
"Otolemur_garnettii", "Pan_paniscus", "Papio_anubis",
"Pongo_abelii"]
for species in val_species:
    Lifted_dfs[species] = pd.read_csv(in_path + species + ".Lifted.tsv", sep='\t',
    low_memory=False, header=None, names=["chr", "start", "end", "label", "score", "strand"])
    Non_lifted_dfs[species] = pd.read_csv(in_path + species + ".Clusters_loc.pep.gz", sep='\t',
    low_memory=False, header=None, names=["chr", "start", "end", "label", "score", "strand"])

#FINAL FILES
final_clusters_info = {}
colnames = [species for species in val_species]
finalnames = ['ClusterID'] + colnames
final_dataset = pd.DataFrame(columns=finalnames)



"""
FUNCTIONS
"""

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def create_dict_with_intervals_for_cluster(list_of_dfs, clusters_dict):
    for species in list_of_dfs:
        for i, row in list_of_dfs[species].iterrows():
            key_str = str(row["label"])
            chr = str(row["chr"])
            start = int(row["start"])
            end = int(row["end"])
            if key_str not in clusters_dict:
                clusters_dict[key_str] = {}
                clusters_dict[key_str][species] = [chr, [start, end]]
            else:
                clusters_dict[key_str][species] = [chr, [start, end]]
    return clusters_dict


def compare_lifted_and_nonlifted(clust_dict1, clust_dict2, final_clust):
    for key in clust_dict1:
        if key in clust_dict1 and key in clust_dict2:
            for species in clust_dict1[key]:
                if species in clust_dict1[key] and species in clust_dict2[key]:
                    list_lifted = clust_dict1[key][species][1]
                    list_nonlifted = clust_dict2[key][species][1]
                    score = getOverlap(list_lifted, list_nonlifted)
                    if key not in final_clust:
                        final_clust[key] = {}
                    if score > 0 and clust_dict1[key][species][0] == clust_dict2[key][species][0]:
                        final_clust[key][species] = 1
                    else:
                        final_clust[key][species] = 0
    return final_clust

def print_final_information(final_clust, final_dataset):
    sorted_keys = sorted(final_clust.keys())
    for key in sorted_keys:
        row = {}
        row["ClusterID"] = key
        for col in final_dataset:
            if col == "ClusterID":
                continue
            elif col in final_clust[key]:
                row[col] = final_clust[key][col]
            else:
                row[col] = np.nan
        final_dataset = final_dataset.append(row, ignore_index = True)
    return final_dataset


#MAIN

clusters_dict_lift = create_dict_with_intervals_for_cluster(Lifted_dfs, clusters_dict_lift)

clusters_dict_nonlift = create_dict_with_intervals_for_cluster(Non_lifted_dfs, clusters_dict_nonlift)

final_clusters_info = compare_lifted_and_nonlifted(clusters_dict_lift,  clusters_dict_nonlift, final_clusters_info)

final_dataset = print_final_information(final_clusters_info, final_dataset)

final_dataset.to_csv(out_file, sep = "\t", index=False)
