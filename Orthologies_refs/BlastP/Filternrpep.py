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
it creating a dictionary where the keys are the IDs without the > symbol and
the values are the nucleotide sequences concatenated
"""


def parse_fasta_create_dictionary_idandsequence_with_nr(fasta_input):
    sequences = {}
    seq = ""
    with gzip.open(fasta_input, "rt") as in_fh:
        for line in in_fh:
            line = line.rstrip()
            if line.startswith(">"):
                if seq != "":
                    if seq in list(sequences.values()):
                        inverted_dict = dict([[v,k] for k,v in sequences.items()])
                        previous_ID = inverted_dict[seq]
                        del sequences[previous_ID]
                        new_ID = previous_ID + "_" + ident
                        sequences[new_ID] = seq
                        seq = ""
                    else:
                        sequences[ident] = seq
                        seq = ""
                ident = line.strip(">")
            else:
                line = line.upper()
                seq += line
        sequences[ident] = seq
    return sequences


"""
In the next function, we are going to print in the output file all those new
proteinIDs after the filtering for each one of the species nr
"""

def write_out_file_for_aasequence_nr(aa_sequences):
    with gzip.open(out_file, "wt") as out_fh:
        for ident in aa_sequences:
                print(">{}\n{}".format(ident, aa_sequences[ident]), file=out_fh)


#MAIN SCRIPT

amino_seqs = parse_fasta_create_dictionary_idandsequence_with_nr(in_file)
write_out_file_for_aasequence_nr(amino_seqs)
