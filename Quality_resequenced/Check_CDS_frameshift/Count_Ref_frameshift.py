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

if len(sys.argv) == 7:
    tree_names = sys.argv[1]
    reference_names = sys.argv[2]
    cds_alignment = sys.argv[3]
    frameshift_file = sys.argv[4]
    gene = sys.argv[5]
    out_file = sys.argv[6]

else:
    sys.exit("The usage shoud be: cds_alignment out_file")



"""
#Print out sorted species file
finalnames = [Species_tag_dict[tag] for tag in Species_tag_dict]
final_cleaned_df = pd.DataFrame(columns=finalnames)
"""

#SPECIES NAMES ASSOCIATES TO TAGS

frameshift_data = pd.read_csv(frameshift_file, sep='\t', low_memory=False)#panda creation
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

refs_dict = {}

alignment_obj = AlignIO.read(cds_alignment, "clustal")
for record in alignment_obj:
    curr = record.id
    if "PD" not in curr:
        refs_dict[curr] = []
    else:
        for key1, value1 in Refs_tag_dict.items():
            if curr.startswith(str(value1)):
                genome = key1
        for key2, value2 in Species_tag_dict.items():
            if value2 == genome:
                refs_dict[key2].append(curr)



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



#Loop to get codons for ALL other primates (using human as reference)
genes_dict = {}
for i, row in frameshift_data.iterrows():
    Gene = row["Gene"]
    if Gene in genes_dict:
        if row["Position"] in genes_dict[Gene]:
            genes_dict[Gene][row["Position"]].append(row['Species'])
        else:
            genes_dict[Gene][row["Position"]] = []
            genes_dict[Gene][row["Position"]].append(row['Species'])
    else:
        genes_dict[Gene] = {}
        genes_dict[Gene][row["Position"]] = []
        genes_dict[Gene][row["Position"]].append(row['Species'])


#print(refs_dict)
#print(genes_dict)


#LOOP TO GET THE perc of affected reseq of a reference

counts_dict = {}
for gene in genes_dict:
    counts_dict[gene] = {}
    for position in genes_dict[gene]:
        counts_dict[gene][position] = {}
        for spc in genes_dict[gene][position]:
            if spc in refs_dict.keys():
                counter = 0
                all_resequenced_mapped = refs_dict[spc]
                if len(all_resequenced_mapped) > 0:
                    for ind in all_resequenced_mapped:
                        if ind in genes_dict[gene][position]:
                            counter += 1
                    perc = counter/len(all_resequenced_mapped)
                    counts_dict[gene][position][spc] = perc
        if not counts_dict[gene][position]:
            del counts_dict[gene][position]

print(counts_dict)
#CREATE PANDAS DATAFRAME

Spc_frame = pd.DataFrame(columns=("Gene", "Position", "Species", "Perc"))
for gene in counts_dict:
    for position in counts_dict[gene]:
        for spc in counts_dict[gene][position]:
            values_to_add = {'Gene': gene, 'Position': round(position, 2), 'Species': spc, "Perc": counts_dict[gene][position][spc]}
            row_to_add = pd.Series(values_to_add)
            Spc_frame = Spc_frame.append(row_to_add, ignore_index = True)
#Spc_frame = pd.DataFrame.from_dict(SpeciesRefs_frameshift_dict, orient='index')
Spc_frame.to_csv(out_file, sep="\t")
        #esto es de prueba
