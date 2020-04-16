#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
from itertools import chain
#import seaborn as sns
from scipy import stats
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import copy
from scipy.stats import sem, t
from scipy import mean
import sys
import re
import os
import gzip


#VARIABLES
#!/usr/bin/env python3

"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""


if len(sys.argv) == 3:
    in_file = sys.argv[1]
    out_file = sys.argv[2]
else:
    sys.exit("The usage shoud be: ./Filternrpep.py in_file output_file")


"""
In this case, we create a function that takes a fasta_file as input, and parse
it creating a dictionary where the keys are the aaseqs without the > symbol and
the values are the IDs (multiple possible IDS for same protein sequence)
"""


def parse_fasta_create_dictionary_idandsequence_with_nr(fasta_input):
    sequences = {}
    seq = ""
    count = 0
    with gzip.open(fasta_input, "rt") as in_fh:
        for line in in_fh:
            count+=1
            print(count)
            line = line.rstrip()
            if line.startswith(">"):
                if seq != "":
                    if seq in sequences.keys():
                        sequences[seq] = sequences[seq] + "_" + ident
                        seq = ""
                    else:
                        sequences[seq] = ident
                        seq = ""
                ident = line.strip(">")
            else:
                line = line.upper()
                seq += line
        if seq in sequences.keys():
            sequences[seq] = sequences[seq] + "_" + ident
        else:
            sequences[seq] = ident
    return sequences


"""
In the next function, we are going to print in the output file all those new
proteinIDs after the filtering for each one of the species nr
"""

def write_out_file_for_aasequence_nr(aa_sequences):
    with gzip.open(out_file, "wt") as out_fh:
        for seq in aa_sequences:
                print(">{}\n{}".format(aa_sequences[seq], seq), file=out_fh)


#MAIN SCRIPT

amino_seqs = parse_fasta_create_dictionary_idandsequence_with_nr(in_file)
write_out_file_for_aasequence_nr(amino_seqs)
