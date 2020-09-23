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

print(sys.argv)

if len(sys.argv) == 6:
    fasta_file = sys.argv[1]
    gene_name = sys.argv[2]
    orthos_file = sys.argv[3]
    human_genes = sys.argv[4]
    out_file = sys.argv[5]
else:
    sys.exit("The usage shoud be: ./Filternrpep.py in_file output_file")


#FUNCTIONS

human_models = pd.read_csv(human_genes, sep='\t', low_memory=False)#panda creation

#DICT SPECIES TAGS
ensembl_ids = human_models['ensembl_peptide_id'].to_list()
gene_names = human_models['hgnc_symbol'].to_list()
ensembl_ids_dict = {ensembl_ids[i]:gene_names[i] for i in range(0, len(ensembl_ids))}

ortho_models = pd.read_csv(orthos_file, sep='\t', low_memory=False)#panda creation

#DICT SPECIES TAGS
species_orthos = ortho_models["Ma's night monkey protein or transcript stable ID"].to_list()
human_orthos = ortho_models["Protein stable ID"].to_list()
orthos_ids_dict = {}
for i in range(0, len(ortho_models)):
    if species_orthos[i] not in orthos_ids_dict:
        orthos_ids_dict[species_orthos[i]] = []
        orthos_ids_dict[species_orthos[i]].append(human_orthos[i])
    else:
        orthos_ids_dict[species_orthos[i]].append(human_orthos[i])



"""
In this case, we create a function that takes a gff_file as input, and parse
it creating a dictionary where the keys are the IDs without the > symbol and
the values are the features of each CDS_id
"""
"""
Reverse complement FUNCTIONS
"""

def reverse(seq):
    """Returns a reversed string"""
    return seq[::-1]


def complement(seq):
    """Returns a complement DNA sequence"""
    complement_dict = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    seq_list = list(seq)
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)


def reverse_complement(seq):
    """"Returns a reverse complement DNA sequence"""
    seq = reverse(seq)
    seq = complement(seq)
    return seq

def select_current_value(current_dict, name):
    for tag in current_dict:
        if name.startswith(tag):
            sel_value = current_dict[tag]
    return sel_value


with gzip.open(fasta_file, "rt") as in_fh, open(out_file, "w") as out_fh :
    for line in in_fh:
        line = line.rstrip()
        if line.startswith(">"):
            ident = line[1:]
            if ident in orthos_ids_dict:
                human_idents = select_current_value(orthos_ids_dict, ident)
                for human_id in human_idents:
                    print(human_id)
                    if str(human_id) in ensembl_ids_dict:
                        gene_curr = select_current_value(ensembl_ids_dict, human_id)
                        if gene_curr == gene_name:
                            print(">"+ident, file=out_fh)
                            new_line = next(in_fh)
                            print(new_line.rstrip(), file=out_fh)
