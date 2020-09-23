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

"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 5:
    gff_file = sys.argv[1]
    orthos_file = sys.argv[2]
    human_genes = sys.argv[3]
    out_file = sys.argv[4]
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
orthos_ids_dict = {species_orthos[i]:human_orthos[i] for i in range(0, len(ensembl_ids))}




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


def parse_gff_file_keep_CDs_ids_in_fasta(gff_input):
    cds = {}
    longest_cds = {}
    strand = {}
    seq = ""
    count = 0
    finalID = ""
    lastID = ""
    with gzip.open(gff_input, "rt") as file1, \
    open(out_file, "wt") as out_fh:
        for line in file1:
            fields1 = line.rstrip().split("\t")
            ProteinID = fields1[8].split(";")[2].split("=")[1].split(",")[0]
            print(">{}\t{}\t{}\t{}\t{}".format(fields1[0], fields1[3],
            fields1[4], fields1[6], ProteinID), file=out_fh)


#MAIN SCRIPT

parse_gff_file_keep_CDs_ids_in_fasta(gff_file)
