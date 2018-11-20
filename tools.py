# miscellaneous tools
import os
import shlex
import subprocess
import sys

import HTSeq
import pysam

PATH = './'

def filter_qual_scores_umi(input_samfile):
    """
    """
    # TODO: filter by quality score - allow only one bp with phred<20 in UMI, and another one in BC.
    af = pysam.AlignmentFile(input_samfile)


def download_genome(genome_dir, ref='hg19'):
    """
    Downloads the hg19 reference genome...
    """
    # TODO: find the hg19 genome???

def preprocess_fastq(fastq1, fastq2, output_dir, **params):
    """
    Performs all the steps before running the alignment. Temporary files
    saved in output_dir.
    """
    
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

    def editd1(seq):
        """Returns the list of sequences with edit distance 1
        It returns a sequence with the same length. So a base insertion
        will result in the last base being truncated"""
        
        seq_len = len(seq)
        
        # Original sequence
        edit_seqs = [seq]
        
        # Insertions
        for i in range(seq_len):
            for b in bases:
                edit_seqs.append(seq[:i] +b + seq[i:-1])
                
        # Deletions
        for i in range(seq_len):
            edit_seqs.append(seq[:i] + seq[i+1:]+'-')
        
        # Nt changes
        for i in range(seq_len):
            for b in bases:
                if b!=seq[i]:
                    edit_seqs.append(seq[:i]+b+seq[i+1:])
        
        # 1nt shift forward
        edit_seqs.append('-'+seq[:-1])
        
        # 1nt shift backward
        edit_seqs.append(seq[1:]+'-')
        
        # Convert Ns to sequences
        output_edit_seqs = []
        for s in edit_seqs:
            output_edit_seqs += convert_degen_seq_to_list(s)
        
        return output_edit_seqs

    def bc_editd1_correction(barcodes,bc_pre='',bc_suf='', bc_start=None, bc_end=None):
        pre_len = len(bc_pre)
        bc_len = len(barcodes[0])
        suf_len = len(bc_suf)
        full_seq_len = pre_len + bc_len + suf_len
        if bc_start is None:
            bc_start = (pre_len-1)
        if bc_end is None:
            bc_end = pre_len + bc_len + 1

        bc_dict = defaultdict(list)
        for bc in barcodes:
            seqs_d1 = unique([s[bc_start:bc_end] for s in unique(editd1(bc_pre + bc + bc_suf))])
            [bc_dict[s].append(bc) for s in seqs_d1]
        bc_map = pd.Series(bc_dict)
        bc_map = bc_map[bc_map.apply(len)==1].apply(lambda s:s[0]).to_dict()
        return bc_map

    def fix_bc(bc,bc_map):
        try:
            bc_fixed = bc_map[bc]
        except:
            bc_fixed = ''
        return bc_fixed

    # Amplicon sequence
    amp_seq = 'NNNNNNNNNNIIIIIIIIGTGGCCGATGTTTCGCATCGGCGTACGACTIIIIIIIIATCCACGTGCTTGAGAGGCCAGAGCATTCGIIIIIIII'

    # Get location of cell barcodes in amplicon:
    bc_len = 8
    bc_starts = []
    c = 0
    while True:
        bc_loc = amp_seq[c:].find('IIIIIIII')
        if bc_loc==-1:
            break
        bc_starts.append(bc_loc + c)
        c = bc_starts[-1] + bc_len

    # Read in barcode sequences
    bc_8nt = pd.read_csv('./barcodes/bc_8nt_v1.csv',names=['barcode'],index_col=0).barcode
    
    # Generate bc_map dictionary for each cell barcode.
    bc3_pre = amp_seq[bc_starts[0]-2:bc_starts[0]]
    bc3_suf = amp_seq[bc_starts[0]+8:bc_starts[0]+10]
    bc3_map = bc_editd1_correction(bc_8nt.values,
                                   bc_pre=bc3_pre,
                                   bc_suf=bc3_suf)

    bc2_pre = amp_seq[bc_starts[1]-2:bc_starts[1]]
    bc2_suf = amp_seq[bc_starts[1]+8:bc_starts[1]+10]
    bc2_map = bc_editd1_correction(bc_8nt.values,
                                   bc_pre=bc2_pre,
                                   bc_suf=bc2_suf)
    bc1_pre = amp_seq[bc_starts[2]-2:bc_starts[2]]
    bc1_suf = 'N'
    bc1_map = bc_editd1_correction(bc_8nt.values,
                                   bc_pre=bc1_pre,
                                   bc_suf=bc1_suf,
                                   bc_start=1,
                                   bc_end=10)

    with gzip.open(fastq1,'rb') as f1, gzip.open(fastq2,'rb') as f2, open(output_dir + 'read1_cellbarcode_headers.fastq','w') as fout:
        while True:
            header2 = f2.readline()
            if len(header2)==0:
                break
            seq2 = f2.readline().decode("utf-8")
            bc1 = fix_bc(seq2[86-1:94],bc1_map)
            bc2 = fix_bc(seq2[48-1:56+1],bc2_map)
            bc3 = fix_bc(seq2[10-1:18+1],bc3_map)
            umi = seq[:10]
            strand2 = f2.readline()
            qual2 = f2.readline()
            
            cellbc_umi = bc3 + bc2 + bc1 +'_' + umi
            header1 = f1.readline().decode("utf-8")
            seq1 = f1.readline().decode("utf-8")
            strand1 = f1.readline().decode("utf-8")
            qual1 = f1.readline().decode("utf-8")
            
            if len(cellbc_umi)==35:
                header1 = '@' + bc3 + bc2 + bc1 +'_' + umi + '_' + qual2.decode("utf-8")[:10] + '_' + header1[1:]
                fout.write(header1)
                fout.write(seq1)
                fout.write(strand1)
                fout.write(qual1)

    return 0

def run_star(genome_dir, output_dir):
    """
    """
    rc = subprocess.call("""STAR --genomeDir {0}/ --runThreadN 8 --readFilesIn {1}/single_cells_barcoded_head.fastq --outFileNamePrefix {1}/single_cells_barcoded_head""".format(genome_dir, output_dir), shell=True)
    return rc

def run_postprocessing(input_dir, output_dir):
    """
    """
