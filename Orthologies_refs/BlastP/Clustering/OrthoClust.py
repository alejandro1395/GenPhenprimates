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


"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 6:
    BBH_path = sys.argv[1]
    selected_species = sys.argv[2]
    tags_file = sys.argv[3]
    out_file = sys.argv[4]
    OMA_out = sys.argv[5]
else:
    sys.exit("The usage shoud be: ./OrthoCLust.py BBH_path tags_file out_file")

#VARIABLES
#!/usr/bin/env python3

#VARIABLES
query_species = ""
ortholog_dict = {}
species_dfs = {}

"""
Read TAGs file foreach primate species used in the analysis
and store it in a Pandas DataFrame
"""

Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames
OrthoClust_df = pd.DataFrame(columns=finalnames)
OMAs = pd.DataFrame(columns=finalnames)

"""
FUNCTIONS
"""

"""
Here, we create a function that concatenate multiple Pandas dataframes
for each species into an allspecies Best Hits dataframe
"""

def concatenate_BBH_all_species_into_dict(path, all_species_dfs):
    dirs = os.listdir(path)
    species = [e for e in dirs if e not in ('qu', 'out')]
    for folder in species:
        path_to_pandas = path + folder + "/" + folder + ".pep.BBH.gz"
        all_species_dfs[folder] = pd.read_csv(path_to_pandas, index_col=None, sep ="\t")
    return all_species_dfs


"""
Then we iterate each row in the BBH concatenated dataframe and keep only the
row if it is the query retrieving the most BBH relationships (whether all
members are in the Cluster reciprocally, send to OMAs)
"""

def create_BBHs_from_rows_in_merged_SBH_data(all_species_dfs, query_species,
bbh_dict, bbh_df, column_names):
    for index, row in all_species_dfs[query_species].iterrows():
        cand_clusters = []
        ortholog_dict[row['Query']] = []
        #retrieve first cluster (corresponding to query row)
        cand_clusters = filter_row_as_candidate_cluster(row, cand_clusters)
        #look for the rest of clusters deriving from it
        for colname in column_names:
            for item in Species_tags['Tag']:
                if str(row[colname]).startswith(item):
                    current_species = Species_tags.loc[Species_tags['Tag'] == item].Species.item()
                    ortholog_dict[row['Query']].append(row[colname])
            target_species = all_species_dfs[current_species].loc[all_species_dfs[current_species]['Query'] == row[colname]]
            cand_clusters = filter_row_as_candidate_cluster(target_species, cand_clusters)
        #IF ORTHO CLUSTER IS an OMA group
        if all(cluster_ind==cand_clusters[0] for cluster_ind in cand_clusters):
            OMAs = append_out_OMAs_pandas_format(ortholog_dict, OMAs, row['Query'])
            OMAs.to_csv(OMA_out, sep ="\t", index=False, mode='a')
        else:

        if not bbh_dict[row['Query']]:
            continue
        else:
            bbh_df = append_out_BBHs_pandas_format(bbh_dict, bbh_df, row['Query'])
    return bbh_df

"""
PROCESS cluster row in order to store ut for comparison with the rest
of candidate reference species
"""


def filter_row_as_candidate_cluster(irow, clusters):
    clust = irow.tolist()
    clust = filter(lambda v: v==v, clust)
    clusters.append(set(clust))
    return clusters

"""
FUnction to parse all pairs of values and print them in the
Query vs multiple targets format
"""

def append_out_row_pandas_format(bbh_dict, bbh_df, key):
    query_row = [key] + bbh_dict[key] + \
    list(np.repeat(np.nan, len(bbh_df.columns)-len(bbh_dict[key])-1))
    bbh_df = bbh_df.append(pd.Series(query_row, index=bbh_df.columns),
    ignore_index=True)
    return bbh_df""



#MAIN

species_dfs = concatenate_SBH_all_species_into_dict(SBH_path, species_dfs)

BBH_df = create_BBHs_from_rows_in_merged_SBH_data(species_dfs, selected_species,
BBH_dict, BBH_df, colnames)

BBH_df.to_csv(out_file, sep = "\t", index=False)
