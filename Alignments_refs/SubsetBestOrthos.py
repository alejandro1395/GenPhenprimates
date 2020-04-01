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
    clusters_path = sys.argv[1]
    tags_file = sys.argv[2]
    out_file = sys.argv[3]
else:
    sys.exit("The usage shoud be: ./OrthoCLust.py clusters_path tags_file out_file")


#VARIABLES

raw_df = pd.read_csv(clusters_path, sep='\t', low_memory=False)#panda creation
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
finalnames = [Species_tag_dict[tag] for tag in Species_tag_dict]
final_cleaned_df = pd.DataFrame(columns=finalnames)




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
FUnction to parse all pairs of values and print them in the
Query vs multiple targets format
"""

def append_out_row_pandas_format(trans_row, output_file):
    trans_row.to_csv(output_file,mode='a',index=False,header=False, sep="\t")

    #Drop nan values for 60-1 species, with clusters missing targets
    #(directed previously to other clusters)

#Reorder based on species names as columns
for index, row in raw_df.iterrows():
    cluster_row = {}
    nans = 0
    for value in row.values.tolist():
        if isNaN(value):
            nans+=1
        else:
            continue
    if nans < 20:
        final_cleaned_df = final_cleaned_df.append(row, ignore_index = True)



final_cleaned_df.to_csv(out_file, sep = "\t", index=False)
