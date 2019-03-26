# make barcode dictionaries
import os
import subprocess
import sys

import pandas as pd
from collections import defaultdict
import gzip
from numpy import unique
import numpy as np
import pickle


#import HTSeq
#import pysam

#PATH = './'
PATH = os.path.dirname(os.path.abspath(__file__))
HOME = os.path.expanduser('~')

print('Path',PATH)

bases = list('ACGT')
def convert_degen_seq_to_list(seq):
    """Uses recursion to convert a degenerate sequence to a list
    For example: AGGN -> [AGGA, AGGC, AGGG, AGGT]"""

    seq_list = []
    N_pos = seq.find('N')
    if N_pos>=0:
        for b in bases:
            seq_list += convert_degen_seq_to_list(seq[:N_pos] + b + seq[N_pos+1:])
    else:
        seq_list.append(seq)
    return seq_list

def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

def generate_edit_dict(barcode_dict):
    edit_dict = {0:defaultdict(list),1:defaultdict(list),2:defaultdict(list),3:defaultdict(list)}
    for bc in barcode_dict:
        for cur_mer in mer8:
            lev_dist = levenshteinDistance(bc,cur_mer)
            if lev_dist==0:
                edit_dict[0][cur_mer].append(bc)
            elif lev_dist==1:
                edit_dict[1][cur_mer].append(bc)
            elif lev_dist==2:
                edit_dict[2][cur_mer].append(bc)
            elif lev_dist==3:
                edit_dict[3][cur_mer].append(bc)
        print(bc,end=' ')
    return edit_dict

mer8 = convert_degen_seq_to_list('NNNNNNNN')

bc_8nt_v1 = pd.read_csv(PATH + '/barcodes/bc_8nt_v1.csv',names=['barcode'],index_col=0).barcode
edit_dict_v1 = generate_edit_dict(bc_8nt_v1)
with open(PATH + '/barcodes/bc_dict_v1.pkl', 'wb') as f:
    pickle.dump(edit_dict_v1, f)

bc_8nt_v2 = pd.read_csv(PATH + '/barcodes/bc_8nt_v2.csv',names=['barcode'],index_col=0).barcode
edit_dict_v2 = generate_edit_dict(bc_8nt_v2)
with open(PATH + '/barcodes/bc_dict_v2.pkl', 'wb') as f:
    pickle.dump(edit_dict_v2, f)

bc_8nt_v3 = pd.read_csv(PATH + '/barcodes/bc_8nt_v3.csv',names=['barcode'],index_col=0).barcode
edit_dict_v3 = generate_edit_dict(bc_8nt_v3)
with open(PATH + '/barcodes/bc_dict_v3.pkl', 'wb') as f:
    pickle.dump(edit_dict_v3, f)