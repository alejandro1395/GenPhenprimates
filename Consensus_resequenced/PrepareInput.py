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
    cds_fasta_path = sys.argv[1]
    tags_file = sys.argv[2]
    out_path = sys.argv[3]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")

#VARIABLES
column1 = "Tag"
column2 = "SpeciesBroad"


#SPECIES NAMES ASSOCIATES TO TAGS
Species_tags = pd.read_csv(tags_file, sep='\t', low_memory=False)#panda creation
colnames = ['Target{}'.format(num) for num in range(1, len(Species_tags))]
finalnames = ['Query'] + colnames

#DICT SPECIES TAGS
individuals_ids = Species_tags['PGDP ID'].to_list()
Species_tags['Species_name'] = Species_tags[['Genus', 'Species']].agg('_'.join, axis=1)
species_names = Species_tags['Species_name'].to_list()
Species_tag_dict = {individuals_ids[i]:species_names[i] for i in range(0, len(individuals_ids))}

"""
#Print out sorted species file
finalnames = [Species_tag_dict[tag] for tag in Species_tag_dict]
final_cleaned_df = pd.DataFrame(columns=finalnames)
"""

"""
FUNCTIONS
"""

"""
Function to check nan values
"""

def isNaN(string):
    return string != string


#Select current species from all species dataframes


def select_current_value(current_dict, name):
    sel_value = [current_dict[tag] for tag in current_dict if name.startswith(tag)][0]
    return sel_value


#Function to select fasta from species file and store it in a dict format

only_genus = False
with open(cds_fasta_path, 'r') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            ident = line[1:]
            species = select_current_value(Species_tag_dict, ident)
            if species.endswith("sp."):
                only_genus = True
            else:
                only_genus = False
                with gzip.open(out_path+species+".cds.fa.gz", 'at') as out_f:
                    print(line, file=out_f)
        elif not only_genus:
            with gzip.open(out_path+species+".cds.fa.gz", 'at') as out_f:
                print(line, file=out_f)
