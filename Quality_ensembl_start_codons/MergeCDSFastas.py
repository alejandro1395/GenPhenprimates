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
    gff_file = sys.argv[2]
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


def parse_gff_file_keep_CDs_ids_in_fasta(gff_input, fasta_input):
    cds = {}
    longest_cds = {}
    strand = {}
    seq = ""
    count = 0
    finalID = ""
    lastID = ""
    with gzip.open(gff_input, "rt") as file1, gzip.open(fasta_input, "rt") as file2:
        for line1, line2 in zip(file1, file2):
            fields1 = line1.rstrip().split("\t")
            curr_strand = fields1[6]
            ProteinID = fields1[8].split(";")[2].split("=")[1].split(",")[0]
            fields2 = line2.rstrip().split("\t")
            strand[ProteinID] = curr_strand
            #CONDITIONS FOR COMPILING CDS SEQUENCES
            if finalID == "":
                finalID = ProteinID
                seq += fields2[1]
            else:
                if ProteinID == finalID:
                    seq += fields2[1]
                else:
                    cds[finalID] = seq
                    finalID = ProteinID
                    seq = ""
                    seq += fields2[1]
            #CONDITION FOR KEEPING ONLY THE LONGEST ONE
            if lastID and ProteinID != lastID:
                longest = ""
                longest = max(cds, key=lambda k: len(cds[k]))
                longest_cds[longest] = cds[longest].upper()
                cds.clear()
            lastID = ProteinID
        cds[finalID] = seq
        longest = ""
        longest = max(cds, key=lambda k: len(cds[k]))
        longest_cds[longest] = cds[longest].upper()
        cds.clear()
    return longest_cds, strand


"""
In the next function, we are going to print in the output file all those new
proteinIDs after the filtering for each one of the species nr
"""

def write_out_file_for_cds(longest, output, strand):
    with gzip.open(output, "wt") as out_fh:
        for ident in longest:
            if strand[ident] == "+":
                print(">{}\n{}".format(ident, longest[ident]), file=out_fh)
            elif strand[ident] == "-":
                reversed = reverse_complement(longest[ident])
                print(">{}\n{}".format(ident, reversed), file=out_fh)



#MAIN SCRIPT

CDS_dict, strand_dict = parse_gff_file_keep_CDs_ids_in_fasta(gff_file, fasta_file)
write_out_file_for_cds(CDS_dict, out_file, strand_dict)
