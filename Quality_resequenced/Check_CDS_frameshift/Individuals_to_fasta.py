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

if len(sys.argv) == 3:
    individual_path = sys.argv[1]
    out_file = sys.argv[2]
else:
    sys.exit("The usage shoud be: ortho_path tags_file pep_fasta_path out_file")



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

individual_name = individual_path.split("\.")[0]
with open(individual_path, 'r') as f:
    for line in f:
        line = line.rstrip()
        fields = line.split(" ")
        ident = fields[0].split("=")[1][:-1]
        print(ident)
        sequence = fields[2]
        with open(out_file, 'a') as out_fh:
            print(">"+ident, file=out_fh)
            print(sequence, file=out_fh)
