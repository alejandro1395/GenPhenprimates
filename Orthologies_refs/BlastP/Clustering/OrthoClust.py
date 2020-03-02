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

if len(sys.argv) == 4:
    BBH_path = sys.argv[1]
    tags_file = sys.argv[2]
    out_file = sys.argv[3]
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

#SPECIES NAMES ASSOCIATES TO TAGS
Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames
#OUTPUT FOLDERS FOR EACH SPECIES
OrthoClust_df = {}
OMAs_df = pd.DataFrame(columns=finalnames)
OMAs_df.to_csv(out_file + "OMAs/all_species.pep.OMAs.gz", sep = "\t", index=False)
target_species = Species_tags['Species'].to_list()
for spc in target_species:
    OrthoClust_df[spc] = pd.DataFrame(columns=finalnames)
    OrthoClust_df[spc].to_csv(out_file + spc + "/" + spc + ".pep.OrthoClust.gz", sep = "\t", index=False)


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
        all_species_dfs[folder] = pd.read_csv(path_to_pandas, index_col=None,
        sep ="\t", low_memory=False)
    return all_species_dfs


"""
Then we iterate each row in the BBH concatenated dataframe and keep only the
row if it is the query retrieving the most BBH relationships (whether all
members are in the Cluster reciprocally, send to OMAs). This functions groups
many other functions used for computing clusters (HOGs-like and OMAs)
"""

def create_clusters_from_rows_in_BBH_data(all_species_dfs,
ortho_df, column_names, omas):
    for query_species in all_species_dfs:
        count = 0
        for index, row in all_species_dfs[query_species].iterrows():
            cand_clusters = {}
            cand_clusters[row['Query']] = []
            count+=1
            print(count)
            #look for the rest of clusters deriving from it
            for colname in column_names:
                if isNaN(row[colname]):
                    continue
                else:
                    curr_species = select_current_species(Species_tags, row[colname])
                    target_species = all_species_dfs[curr_species].loc[all_species_dfs[curr_species]['Query'] == row[colname]]
                    cand_clusters = filter_row_as_candidate_cluster(target_species, cand_clusters)
                cand_clusters[row['Query']].append(row[colname])
            #IF ORTHO CLUSTER IS an OMA group
            cand_clusters[row['Query']] = set(cand_clusters[row['Query']])
            if all(set([key] + list(cand_clusters[key]))==set([row['Query']] + list(cand_clusters[row['Query']])) for key in cand_clusters):
                append_out_row_pandas_format(cand_clusters, omas, row['Query'], out_file + "OMAs/all_species.pep.OMAs.gz")
                all_species_dfs = set_nan_values_for_included_orthologous_proteins(all_species_dfs,
                [row['Query']], cand_clusters[row['Query']])
            #whether not check the longest one in terms of retrieval of orthologies
            else:
                max_length, max_key = GetMaxFlox(cand_clusters)
                curr_species = select_current_species(Species_tags, max_key)
                append_out_row_pandas_format(cand_clusters,
                ortho_df[curr_species], max_key, out_file + curr_species + "/" + curr_species + ".pep.OrthoClust.gz")
                #print(ortho_df[curr_species])
                all_species_dfs = set_nan_values_for_included_orthologous_proteins(all_species_dfs,
                [row['Query']], cand_clusters[row['Query']])
    return ortho_df, omas


"""
Function to check nan values
"""

def isNaN(string):
    return string != string


"""
PROCESS cluster row in order to store it for comparison with the rest
of candidate reference species
"""


def filter_row_as_candidate_cluster(irow, clusters):
    clust = irow.values.tolist()[0][1:]
    clust = filter(lambda v: v==v, clust)
    clusters[str(irow.values.tolist()[0][0])] = set(clust)
    return clusters

"""
FUnction to parse all pairs of values and print them in the
Query vs multiple targets format
"""

def append_out_row_pandas_format(bbh_dict, bbh_df, key, output_file):
    query_row = [key] + list(bbh_dict[key]) + \
    list(np.repeat(np.nan, len(bbh_df.columns)-len(bbh_dict[key])-1))
    new_row = pd.Series(query_row, index=bbh_df.columns)
    trans_row = new_row.to_frame().T
    trans_row.to_csv(output_file,mode='a',index=False,header=False, sep="\t")

"""
Retrieve longest dict element
"""

def GetMaxFlox(flows):
    maks=max(flows, key=lambda k: len(flows[k]))
    return len(flows[maks]), maks

"""
Select current species from all species dataframes
"""

def select_current_species(species_tags, name):
    for item in species_tags['Tag']:
        if str(name).startswith(str(item)):
            current_species = species_tags.loc[species_tags['Tag'] == item].Species.item()
    return current_species


"""
Replace entries in all dataframes for which the ID of the protein has already
been included in almost one ortholog cluster group before
"""
def set_nan_values_for_included_orthologous_proteins(df, query, query_values):
    val_to_replace = query + list(query_values)
    for spc in df:
        df[spc] = pd.DataFrame(np.where(df in val_to_replace, np.nan, df[spc]), index=df[spc].index, columns=df[spc].columns)
    return df


"""
Print output for each one of the dataframes storing
orthologous groups
"""

def print_final_output_for_OMAs_and_ref_groups(omas, orthos_df, out_file):
    omas.to_csv(out_file + "OMAs/all_species.pep.OMAs.gz", sep = "\t", index=False)
    for spc in orthos_df:
        orthos_df[spc].to_csv(out_file + spc + "/" + spc + ".pep.OrthoClust.gz", sep = "\t", index=False)


#MAIN

species_dfs = concatenate_BBH_all_species_into_dict(BBH_path, species_dfs)

OrthoClust_df, OMAs_df = create_clusters_from_rows_in_BBH_data(species_dfs,
OrthoClust_df, colnames, OMAs_df)

#print_final_output_for_OMAs_and_ref_groups(OMAs_df, OrthoClust_df, out_file)
