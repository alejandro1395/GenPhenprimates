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
    SBH_path = sys.argv[1]
    selected_species = sys.argv[2]
    tags_file = sys.argv[3]
    out_file = sys.argv[4]
else:
    sys.exit("The usage shoud be: ./FinderBBH.py SBH_path tags_file out_file")

#VARIABLES
#!/usr/bin/env python3

#VARIABLES
query_species = ""
BBH_dict = {}
species_dfs = {}

"""
Read TAGs file foreach primate species used in the analysis
and store it in a Pandas DataFrame
"""

Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames
BBH_df = pd.DataFrame(columns=finalnames)

"""
FUNCTIONS
"""

"""
Here, we create a function that concatenate multiple Pandas dataframes
for each species into an allspecies Best Hits dataframe
"""

def concatenate_SBH_all_species_into_dict(path, all_species_dfs):
    files = os.walk(path)
    for file_in in files:
        for spc in file_in[2]:
            print(spc)
            path_to_pandas = path + spc
            all_species_dfs[spc.split(".")[0]] = pd.read_csv(path_to_pandas, index_col=None, sep ="\t")
    return all_species_dfs


"""
Then we iterate each row in the SBH concatenated dataframe and keep only the
values in each row that are reciprocal best hits
"""

def create_BBHs_from_rows_in_merged_SBH_data(all_species_dfs, query_species,
bbh_dict, bbh_df, column_names):
    for index, row in all_species_dfs[query_species].iterrows():
        bbh_dict[row['Query']] = []
        for colname in column_names:
            if isNaN(row[colname]):
                continue
            else:
                for item in Species_tags['Tag']:
                    if str(row[colname]).startswith(item):
                        current_species = Species_tags.loc[Species_tags['Tag'] == item].Species.item()
                target_species = all_species_dfs[current_species].loc[all_species_dfs[current_species]['Query'] == row[colname]]
                if row['Query'] in target_species.values:
                    bbh_dict[row['Query']].append(row[colname])
        if not bbh_dict[row['Query']]:
            continue
        else:
            bbh_df = append_out_BBHs_pandas_format(bbh_dict, bbh_df, row['Query'])
    return bbh_df

"""
Function to check nan values
"""

def isNaN(string):
    return string != string

"""
FUnction to parse all BBH pairs of values and print them in the
Query vs multiple targets format
"""

def append_out_BBHs_pandas_format(bbh_dict, bbh_df, key):
    query_row = [key] + bbh_dict[key] + \
    list(np.repeat(np.nan, len(bbh_df.columns)-len(bbh_dict[key])-1))
    bbh_df = bbh_df.append(pd.Series(query_row, index=bbh_df.columns),
    ignore_index=True)
    return bbh_df

#MAIN

species_dfs = concatenate_SBH_all_species_into_dict(SBH_path, species_dfs)

BBH_df = create_BBHs_from_rows_in_merged_SBH_data(species_dfs, selected_species,
BBH_dict, BBH_df, colnames)

BBH_df.to_csv(out_file, sep = "\t", index=False)
