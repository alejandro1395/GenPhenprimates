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
    positions_file = sys.argv[2]
    gene_name = sys.argv[3]
    out_file = sys.argv[4]

else:
    sys.exit("The usage shoud be: cds_alignment out_file")



"""
FUNCTIONS
"""

stop_codons = ["TGA", "TAG", "TAA"]
stop_codon_info = pd.read_csv(positions_file, sep='\t', low_memory=False,
names=["Index", "Gene", "Species", "Position", "Perc_Reseq"])#panda creation


#REFS_DICT AND RESEQUENCED_DICT
refs_dict = {}
refs_pos = {}
alignment_obj = AlignIO.read(cds_alignment, "clustal")
for record in alignment_obj:
    curr = str(record.id)
    seq = str(record.seq).upper()
    refs_dict[curr] = seq
    positions = [m.start(0) for m in re.finditer("(A|G|C|T)", str(record.seq), flags=re.IGNORECASE)]
    refs_pos[curr] = positions


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

#Create dict with positions
genes_dict = {}
for i, row in stop_codon_info.iterrows():
    Species = row["Species"]
    if Species in genes_dict:
        if row['Position'] in genes_dict[Species]:
            genes_dict[Species][row['Position']].append(row['Perc_Reseq'])
        else:
            genes_dict[Species][row['Position']] = []
            genes_dict[Species][row['Position']].append(row['Perc_Reseq'])
    else:
        genes_dict[Species] = {}
        genes_dict[Species][row['Position']] = []
        genes_dict[Species][row['Position']].append(row['Perc_Reseq'])



#Loop through cases taking alignment of references into consideration

Spc_frame = pd.DataFrame(columns=("Gene", "Species", "Position", "Perc_Reseq",
"Perc_Refs", "Pos_rel"))

for spc in genes_dict:
    for pos in genes_dict[spc]:
        curr_positions = refs_pos[spc]
        curr_sequence = refs_dict[spc]
        pos_rel = pos/len(curr_positions)
        curr_codon = get_codon_info(curr_sequence, pos-1, curr_positions)
        alignment_position = curr_positions[pos-1]
        ref_codons_dict = {}
        for ref in refs_dict:
            if ref != spc:
                ref_codon = get_codon_info(refs_dict[ref], pos-1, curr_positions)
                ref_codons_dict[ref] = ref_codon
        count_validated = 0
        for ref_spc in ref_codons_dict:
            if ref_codons_dict[ref_spc] == curr_codon:
                count_validated += 1
        perc_refs = count_validated/(len(refs_dict)-1)
        values_to_add = {'Gene': gene_name, 'Species': spc,
        'Position': pos, "Perc_Reseq": round(genes_dict[spc][pos][0],2),
        "Perc_Refs": round(perc_refs, 2), "Pos_rel": round(pos_rel, 2)}
        row_to_add = pd.Series(values_to_add)
        Spc_frame = Spc_frame.append(row_to_add, ignore_index = True)




"""
stop_codon_pos_ref = {}
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
                    stop_codon_pos_ref[index] = 0
                    for spc in reseq_codons_dict:
                        if reseq_codons_dict[spc] in ["TGA", "TAG", "TAA"]:
                            stop_codon_pos_ref[index] += 1
                    stop_codon_pos_ref[index] = stop_codon_pos_ref[index]/(len(reseq_codons_dict))
"""


if "Homo_sapiens" in refs_dict:
    Spc_frame.to_csv(out_file, sep="\t")
