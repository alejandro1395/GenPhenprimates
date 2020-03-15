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

if len(sys.argv) == 5:
    BBH_path = sys.argv[1]
    tags_file = sys.argv[2]
    traits_file = sys.argv[3]
    out_file = sys.argv[4]
else:
    sys.exit("The usage shoud be: ./OrthoCLust.py BBH_path tags_file traits_file out_file")

#VARIABLES
#!/usr/bin/env python3

#VARIABLES
query_species = ""
ortholog_dict = {}
species_dfs = []
used_genes = set()

"""
Read TAGs file foreach primate species used in the analysis
and store it in a Pandas DataFrame
"""

#SPECIES NAMES ASSOCIATES TO TAGS
Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames

#TRAIT COUTNS ASSOCIATED TO SPECIES
Trait_covs = pd.read_csv(traits_file, sep='\t', low_memory=False)#panda creation

#OUTPUT FOLDERS FOR EACH SPECIES
OrthoClust_df = {}
OMAs_df = pd.DataFrame(columns=finalnames)
OMAs_df.to_csv(out_file + "OMAs/all_species.pep.OMAs.gz", sep = "\t", index=False)
target_species = Species_tags['Species'].to_list()
for spc in target_species:
    OrthoClust_df[spc] = pd.DataFrame(columns=finalnames)
    OrthoClust_df[spc].to_csv(out_file + spc + "/" + spc + ".pep.OrthoClustHumanprior_newBBH.gz", sep = "\t", index=False)


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
        df_species = pd.read_csv(path_to_pandas, index_col=None,
        sep ="\t", low_memory=False)
        all_species_dfs.append(df_species)
    return all_species_dfs

def concatenate_into_single_df(all_species_dfs):
    result = pd.concat(all_species_dfs)
    return result

"""
Then we iterate each row in the BBH concatenated dataframe and keep only the
row if it is the query retrieving the most BBH relationships (whether all
members are in the Cluster reciprocally, send to OMAs). This functions groups
many other functions used for computing clusters (HOGs-like and OMAs)
"""

def create_clusters_from_rows_in_BBH_data(all_species_dfs,
ortho_df, column_names, omas, used_genes):
    count = 0
    for index, row in all_species_dfs.iterrows():
        #candidate clusters to select the longest one
        cand_clusters = {}
        #current genes to select the ones not to take into account (STRICT ANALYSIS AGAINST DUPLICATIONS)
        current_genes = set()
        if isNaN(row['Query']):
            continue
        else:
            cand_clusters[row['Query']] = []
            current_genes.add(row['Query'])
            count+=1
            print(count)
            #look for the rest of clusters deriving from it
            for colname in column_names:
                if isNaN(row[colname]):
                    continue
                else:
                    #curr_species = select_current_species(Species_tags, row[colname])
                    target_species = all_species_dfs.loc[all_species_dfs['Query'] == row[colname]]
                    cand_clusters, current_genes = filter_row_as_candidate_cluster(target_species, cand_clusters,
                    current_genes)
                current_genes.add(row[colname])
                cand_clusters[row['Query']].append(row[colname])
            #IF ORTHO CLUSTER IS an OMA group
            cand_clusters[row['Query']] = set(cand_clusters[row['Query']])
            if any(val in list(used_genes) for val in list(current_genes)):
                continue
            else:
                max_length, max_key, curr_species = GetMaxFlox(cand_clusters, Trait_covs, Species_tags)
                used_list = [max_key] + list(cand_clusters[max_key])
                append_out_row_pandas_format(cand_clusters,
                ortho_df[curr_species], max_key, out_file + curr_species + "/" + curr_species + ".pep.OrthoClustHumanprior_newBBH.gz")
            for value in list(current_genes):
                used_genes.add(value)
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

def filter_row_as_candidate_cluster(irow, clusters, curr_genes):
    clust = irow.values.tolist()[0][1:]
    clust = [x for x in clust if str(x) != 'nan']
    clusters[str(irow.values.tolist()[0][0])] = set(clust)
    curr_genes.add(str(irow.values.tolist()[0][0]))
    for element in clust:
        curr_genes.add(element)
    return clusters, curr_genes

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
Select current species from all species dataframes
"""

def select_current_species(species_tags, name, col):
    print(name)
    sel_species = species_tags.loc[species_tags[col].str.startswith(str(name), na=False)].Species.item()
    print(sel_species)
    return sel_species


"""
Retrieve longest dict element
"""

def GetMaxFlox(l, Trait_covs, species_tag):
    maks=max(l, key=lambda k: len(l[k]))
    maxlist = [s for s in l if len(l[s]) == len(l[maks])]
    #added filter to check human presence
    if any(key.startswith("ENSP00") for key in maxlist):
        highest_spcs = [i for i in maxlist if i.startswith("ENSP00")][0]
        final_spc = select_current_species(species_tag, highest_spcs)
    #when no human, priotitize species with much more trait information
    else:
        highest_num = 0
        highest_spcs = ""
        for i in maxlist:
            column = "Tag"
            curr = select_current_species(species_tag, i, column)
            column = "SpeciesBroad"
            num_traits = select_current_species(Trait_covs, curr, column)
            if (num_traits > highest_num):
                highest_num = num_traits
                highest_spcs = i
        final_spc = select_current_species(species_tag, highest_spcs)
    return len(l[highest_spcs]), highest_spcs, final_spc



"""
Replace entries in all dataframes for which the ID of the protein has already
been included in almost one ortholog cluster group before
"""
def set_nan_values_for_other_ref_orthologous_proteins(df, query, query_values, species_tags):
    val_to_replace = query + list(query_values)
    df.values[df.applymap(str).isin(val_to_replace)] = np.nan



"""
Print output for each one of the dataframes storing
orthologous groups
"""

def print_final_output_for_OMAs_and_ref_groups(omas, orthos_df, out_file):
    omas.to_csv(out_file + "OMAs/all_species.pep.OMAs.gz", sep = "\t", index=False)
    for spc in orthos_df:
        orthos_df[spc].to_csv(out_file + spc + "/" + spc + ".pep.OrthoClust_newBBH_int.gz", sep = "\t", index=False)



#######
#MAIN##
#######

species_dfs = concatenate_BBH_all_species_into_dict(BBH_path, species_dfs)

concatenated_df = concatenate_into_single_df(species_dfs)

OrthoClust_df, OMAs_df = create_clusters_from_rows_in_BBH_data(concatenated_df,
OrthoClust_df, colnames, OMAs_df, used_genes)

#print_final_output_for_OMAs_and_ref_groups(OMAs_df, OrthoClust_df, out_file)
