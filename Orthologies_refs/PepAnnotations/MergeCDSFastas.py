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


if len(sys.argv) == 4:
    fasta_file = sys.argv[1]
    gff_file =sys.argv[2]
    out_file = sys.argv[3]
else:
    sys.exit("The usage shoud be: ./Filternrpep.py in_file output_file")


#FUNCTIONS

"""
In this case, we create a function that takes a gff_file as input, and parse
it creating a dictionary where the keys are the IDs without the > symbol and
the values are the features of each CDS_id
"""


def parse_gff_file_keep_CDs_ids_in_fasta(gff_input, fasta_input):
    cds = {}
    seq = ""
    count = 0
    finalID = ""
    with gzip.open(gff_input, "rt") as file1, gzip.open(fasta_input, "rt") as file2:
        for line1, line2 in zip(file1, file2):
            fields1 = line1.rstrip().split("\t")
            ID = fields1[8].split(";")[0].split("=")[1]
            fields2 = line2.rstrip().split("\t")
            if finalID == "":
                finalID = ID
                seq += fields2[1]
            else:
                if ID == finalID:
                    seq += fields2[1]
                else:
                    cds[finalID] = seq
                    finalID = ID
                    seq = ""
                    seq += fields2[1]
        cds[finalID] = seq
    return cds


"""
In the next function, we are going to print in the output file all those new
proteinIDs after the filtering for each one of the species nr
"""

def write_out_file_for_cds(cds):
    with gzip.open(out_file, "wt") as out_fh:
        for ident in cds:
            print(">{}\n{}".format(ident, cds[ident]), file=out_fh)




#MAIN SCRIPT

CDS_dict = parse_gff_file_keep_CDs_ids_in_fasta(gff_file, fasta_file)
write_out_file_for_cds(CDS_dict)
