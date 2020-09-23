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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import argparse
import glob



"""
Here we create the function for checking the input parameters and saving
them in different variables, error if the usage is not good
"""

if len(sys.argv) == 5:
    cds_alignment = sys.argv[1]
    gene_name = sys.argv[2]
    spc_name = sys.argv[3]
    out_file = sys.argv[4]

else:
    sys.exit("The usage shoud be: cds_alignment out_file")



"""
FUNCTIONS
"""

stop_codons = ["TGA", "TAG", "TAA"]

#REFS_DICT AND RESEQUENCED_DICT
refs_dict = {}
resequenced_dict = {}
refs_pos = {}
resequenced_pos = {}
alignment_obj = AlignIO.read(cds_alignment, "clustal")
for record in alignment_obj:
    curr = record.id
    seq = str(record.seq)
    if "PD" not in curr:
        refs_dict[curr] = seq
        positions = [m.start(0) for m in re.finditer("(A|G|C|T)", str(record.seq), flags=re.IGNORECASE)]
        refs_pos[curr] = positions
    else:
        resequenced_dict[curr] = seq
        positions = [m.start(0) for m in re.finditer("(A|G|C|T)", str(record.seq), flags=re.IGNORECASE)]
        resequenced_pos[curr] = positions


#first_pos = re.search("(A|G|C|T)", str(human_aligned_cds), flags=re.IGNORECASE).start()
#last_pos =  re.search("(A|G|C|T)", str(human_aligned_cds)[-1], flags=re.IGNORECASE).end()




"""
get codon information
"""

def get_codon_info(sequence_aligned, index_pos, positions):
    codon_seq = sequence_aligned[positions[index_pos]] + sequence_aligned[positions[index_pos+1]] + \
    sequence_aligned[positions[index_pos+2]]
    return codon_seq


"""
remove empty keys from positions python
"""

def remove_empty_keys(d):
    for k in d.keys():
        if not d[k]:
            del d[k]





"""
Select current species from all species dataframes
"""

def select_current_value(current_dict, name):
    sel_value = [current_dict[tag] for tag in current_dict if name.startswith(str(tag))][0]
    return sel_value


#MAIN

#Search for STOP codons

stop_codon_pos_ref = {}
stop_codon_pos_frame= {}
for reseq in resequenced_pos:
    positions = resequenced_pos[reseq]
    for i in range(0, len(positions), 3):
        curr_seq = resequenced_dict[reseq]
        curr_seq = curr_seq.upper()
        if i < len(positions)-3:
            curr_codon = get_codon_info(curr_seq, i, positions)
            if curr_codon in ["TGA", "TAG", "TAA"]:
                reseq_codons = []
                reseq_codons_dict = {}
                for record in alignment_obj:
                    seq_species = str(record.seq)
                    seq_species = seq_species.upper()
                    if record.id in refs_dict:
                        ref_codon = get_codon_info(seq_species, i, positions)
                        ref_positions = refs_pos[record.id]
                    else:
                        reseq_codon = get_codon_info(seq_species, i, positions)
                        reseq_codons.append(reseq_codon)
                        reseq_codons_dict[record.id] = reseq_codon
                #print(miss)
                if any(cod in stop_codons for cod in reseq_codons) and not i == len(positions)-3 \
                and ref_codon not in stop_codons and "-" not in ref_codon:
                    index_alignment = positions[i]
                    index = int(ref_positions.index(index_alignment))+1
                    if int(ref_positions.index(index_alignment)) % 3 == 0:
                        stop_codon_pos_frame[index] = "Yes"
                    else:
                        stop_codon_pos_frame[index] = "No"
                    stop_codon_pos_ref[index] = 0
                    for spc in reseq_codons_dict:
                        if reseq_codons_dict[spc] in ["TGA", "TAG", "TAA"]:
                            stop_codon_pos_ref[index] += 1
                    stop_codon_pos_ref[index] = stop_codon_pos_ref[index]/(len(reseq_codons_dict))


#CREATE PANDAS DATAFRAME

Spc_frame = pd.DataFrame(columns=("Gene", "Species", "Position", "Perc", "Frameshift"))
for pos in stop_codon_pos_ref:
    values_to_add = {'Gene': gene_name, 'Species': spc_name,
    'Position': pos, "Perc": stop_codon_pos_ref[pos], "Frameshift": stop_codon_pos_frame[pos]}
    row_to_add = pd.Series(values_to_add)
    Spc_frame = Spc_frame.append(row_to_add, ignore_index = True)
#Spc_frame = pd.DataFrame.from_dict(SpeciesRefs_frameshift_dict, orient='index')

Spc_frame.to_csv(out_file, sep="\t")
