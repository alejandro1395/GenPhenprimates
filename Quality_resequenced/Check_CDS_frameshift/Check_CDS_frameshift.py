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

if len(sys.argv) == 6:
    tree_names = sys.argv[1]
    reference_names = sys.argv[2]
    cds_alignment = sys.argv[3]
    gene = sys.argv[4]
    out_file = sys.argv[5]

else:
    sys.exit("The usage shoud be: cds_alignment out_file")



"""
#Print out sorted species file
finalnames = [Species_tag_dict[tag] for tag in Species_tag_dict]
final_cleaned_df = pd.DataFrame(columns=finalnames)
"""

#SPECIES NAMES ASSOCIATES TO TAGS

Species_tags = pd.read_csv(tree_names, sep='\t', low_memory=False)#panda creation
Ref_tags = pd.read_csv(reference_names, sep='\t', low_memory=False)

#DICT SPECIES TAGS
Genomes_names = Species_tags['Genomes_names'].to_list()
Tree_names = Species_tags['Tree_names'].to_list()
Species_tag_dict = {Tree_names[i]:Genomes_names[i] for i in range(0, len(Genomes_names))}
Reference_species = Ref_tags['Species'].to_list()
Reference_idents = Ref_tags['Tag'].to_list()
Refs_tag_dict = {Reference_species[i]:Reference_idents[i] for i in range(0, len(Reference_species))}


"""
FUNCTIONS
"""

alignment_obj = AlignIO.read(cds_alignment, "clustal")
species_alignment = [record.id for record in alignment_obj]


#Let's iterate and take the human sequence as our reference
human_aligned_cds = [record.seq for record in alignment_obj if record.id == "Homo_sapiens"][0]

#print(human_aligned_cds[3990:3993])
#Then we search the first index of its sequence
#From that index we start to calculate its open-reading-frame

#first_pos = re.search("(A|G|C|T)", str(human_aligned_cds), flags=re.IGNORECASE).start()
#last_pos =  re.search("(A|G|C|T)", str(human_aligned_cds)[-1], flags=re.IGNORECASE).end()

positions = [m.start(0) for m in re.finditer("(A|G|C|T)", str(human_aligned_cds), flags=re.IGNORECASE)]
first_pos = positions[0]
last_pos = positions[-1]



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


def check_whether_ref_with_resequenced(record_info, species_alignment):
    curr_spc = record_info.id
    if curr_spc in Tree_names:
        case = "Ref"
        genome_info = select_current_value(Species_tag_dict, record_info.id)
        ident_info = select_current_value(Refs_tag_dict, genome_info)
        if any(spc.startswith(str(ident_info)) for spc in species_alignment) and not "PD" in curr_spc:
            include = True
        else:
            include = False
    else:
        case = "Reseq"
        tree_info = None
        for key, value in Refs_tag_dict.items():
            if curr_spc.startswith(str(value)):
                genome_info = key
        print(genome_info)
        for key, value in Species_tag_dict.items():
            if value == genome_info:
                tree_info = key
        if tree_info is not None:
            if tree_info in species_alignment:
                include = True
            else:
                include = False
        else:
            include = False
    return case, include




"""
Select current species from all species dataframes
"""

def select_current_value(current_dict, name):
    sel_value = [current_dict[tag] for tag in current_dict if name.startswith(str(tag))][0]
    return sel_value





#loop to get codons for ALL other primates (using human as reference)
SpeciesRefs_frameshift_dict = {}
for i in range(0, len(positions), 3):
    curr_codons = []
    cases_dict = {}
    curr_codons_dict = {}
    previos_codons_dict = {}
    next_codons_dict = {}
    #Now I have to look the rpevious codon and the final codon , appart from
    #the affected codon
    for record in alignment_obj:
        case_spc, include_spc = check_whether_ref_with_resequenced(record, species_alignment)
        if include_spc:
            seq_species = record.seq
            codon = get_codon_info(seq_species, i, positions)
            curr_codons.append(codon)
            curr_codons_dict[record.id] = codon
            cases_dict[record.id] = case_spc
            #previous
            jump = i-3
            if jump < 3:
                previos_codons_dict[record.id] = "---"
            else:
                while jump >= 3:
                    jump -= 3
                    print(jump)
                    previous_codon = get_codon_info(seq_species, jump, positions)
                    if previous_codon != "---":
                        break
                previos_codons_dict[record.id] = previous_codon
            #next
            jump = i+3
            if jump > len(positions)-3:
                next_codons_dict[record.id] = "---"
            else:
                while jump < len(positions)-3:
                    jump += 3
                    next_codon = get_codon_info(seq_species, jump, positions)
                    if next_codon != "---":
                        break
                next_codons_dict[record.id] = next_codon
    miss = 0
    for element in curr_codons:
        if str(element) == "---":
            miss += 1
    #print(miss)
    if miss/len(curr_codons) > 0:
        continue
    else:
        if any("-" in str(cod) for cod in curr_codons):
            perc_frame = (i+1)/(len(positions)+1)*100
            SpeciesRefs_frameshift_dict[perc_frame] = []
            for spc in curr_codons_dict:
                if "-" in str(curr_codons_dict[spc]) and not "-" in previos_codons_dict[spc] \
                and not "-" in next_codons_dict[spc]:
                    SpeciesRefs_frameshift_dict[perc_frame].append(spc)


#CREATE PANDAS DATAFRAME
SpeciesRefs_frameshift_dict_final = {key: value for key, value in SpeciesRefs_frameshift_dict.items() if value}

Spc_frame = pd.DataFrame(columns=("Gene", "Position", "Species"))
for pos in SpeciesRefs_frameshift_dict_final:
    for spc in SpeciesRefs_frameshift_dict_final[pos]:
        values_to_add = {'Gene': gene, 'Position': round(pos, 2), 'Species': spc}
        row_to_add = pd.Series(values_to_add)
        Spc_frame = Spc_frame.append(row_to_add, ignore_index = True)
#Spc_frame = pd.DataFrame.from_dict(SpeciesRefs_frameshift_dict, orient='index')
Spc_frame.to_csv(out_file, sep="\t")
        #esto es de prueba
