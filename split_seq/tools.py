# miscellaneous tools
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
PATH = os.path.dirname(__file__)
HOME = os.path.expanduser('~')

STAR_PATH = os.path.join(HOME, 'split_seq_reqs', 'bin', 'STAR')
if not os.path.exists(STAR_PATH):
    STAR_PATH = 'STAR'

SAMTOOLS_PATH = os.path.join(HOME, 'split_seq_reqs', 'bin', 'samtools')
if not os.path.exists(SAMTOOLS_PATH):
    SAMTOOLS_PATH = 'samtools'



def download_genome(genome_dir, ref='hg19'):
    """
    Downloads the hg19 reference genome...
    """
    # TODO: find the hg19 genome???
    
def make_combined_genome(species, fasta_filenames, output_dir):
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create a combined fasta file with species names added to the start of each chromosome name
    cur_fa = fasta_filenames[0]
    cur_species = species[0]
    if fasta_filenames[0].split('.')[-1]=='gz':
        command = """gunzip -cd {0} | awk 'substr($0,1,1)==">"{{print ">{1}_"substr($1,2,length($1)-1),$2,$3,$4}}substr($0,1,1)!=">"{{print $0}}' > {2}/genome.fa""".format(cur_fa, cur_species, output_dir)
    else:
        command = """cat {0} | awk 'substr($0,1,1)==">"{{print ">{1}_"substr($1,2,length($1)-1),$2,$3,$4}}substr($0,1,1)!=">"{{print $0}}' > {2}/genome.fa""".format(cur_fa, cur_species, output_dir)
    rc = subprocess.call(command, shell=True)
    
    for i in range(1,len(species)):
        cur_fa = fasta_filenames[i]
        cur_species = species[i]
        if fasta_filenames[0].split('.')[-1]=='gz':
            command = """gunzip -cd {0} | awk 'substr($0,1,1)==">"{{print ">{1}_"substr($1,2,length($1)-1),$2,$3,$4}}substr($0,1,1)!=">"{{print $0}}' >> {2}/genome.fa""".format(cur_fa, cur_species, output_dir)
        else:
            command = """cat {0} | awk 'substr($0,1,1)==">"{{print ">{1}_"substr($1,2,length($1)-1),$2,$3,$4}}substr($0,1,1)!=">"{{print $0}}' >> {2}/genome.fa""".format(cur_fa, cur_species, output_dir)
        rc = subprocess.call(command, shell=True)
        
def split_attributes(s):
    """ Returns a dictionary from string of attributes in a GTF/GFF file
    """
    att_list = s[:-1].split('; ')
    att_keys = [a.split(' ')[0] for a in att_list]
    att_values = [' '.join(a.split(' ')[1:]) for a in att_list]
    return dict(zip(att_keys,att_values))

def get_attribute(s,att):
    att_value = ''
    try:
        att_value = split_attributes(s)[att].strip('"')
    except:
        att_value = ''
    return att_value
        
def make_gtf_annotations(species, gtf_filenames, output_dir, splicing):
    splicing = splicing=='True'

    # Load the GTFs
    names = ['Chromosome',
         'Source',
         'Feature',
         'Start',
         'End',
         'Score',
         'Strand',
         'Frame',
         'Attributes']

    gtfs = {}
    for i in range(len(species)):
        s = species[i]
        filename = gtf_filenames[i]
        gtfs[s] = pd.read_csv(filename,sep='\t',names=names,comment='#',engine='python')
    
    # TODO: allow users to specify the gene biotypes that they want to keep
    # For now we keep the following
    gene_biotypes_to_keep = ['protein_coding',
                             'lincRNA',
                             'antisense',
                             'IG_C_gene',
                             'IG_C_pseudogene',
                             'IG_D_gene',
                             'IG_J_gene',
                             'IG_J_pseudogene',
                             'IG_V_gene',
                             'IG_V_pseudogene',
                             'TR_C_gene',
                             'TR_D_gene',
                             'TR_J_gene',
                             'TR_J_pseudogene',
                             'TR_V_gene',
                             'TR_V_pseudogene']
    if splicing:
        # Generate a combined GTF with only the gene annotations
        gtf_gene_combined = gtfs[species[0]].query('Feature=="gene"')
        gtf_gene_combined.loc[:,'Chromosome'] = species[0] + '_' + gtf_gene_combined.Chromosome.apply(lambda s:str(s))
        for i in range(1,len(species)):
            gtf_gene_combined_temp = gtfs[species[i]].query('Feature=="gene"')
            gtf_gene_combined_temp.loc[:,'Chromosome'] = species[i] + '_' + gtf_gene_combined_temp.Chromosome.apply(lambda s:str(s))
            gtf_gene_combined = pd.concat([gtf_gene_combined,gtf_gene_combined_temp])
        gene_biotypes = gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_biotype'))
        #gtf_gene_combined = gtf_gene_combined.iloc[np.where(gene_biotypes.isin(gene_biotypes_to_keep).values)]
        gtf_gene_combined.index = range(len(gtf_gene_combined))
        gtf_gene_combined.to_csv(output_dir + '/genes.gtf',sep='\t',index=False)
    
    # Generate a combined GTF with only the exon annotations
    gtf_exon_combined = gtfs[species[0]].query('Feature=="exon"')
    gtf_exon_combined.loc[:,'Chromosome'] = species[0] + '_' + gtf_exon_combined.Chromosome.apply(lambda s:str(s))
    for i in range(1,len(species)):
        gtf_exon_combined_temp = gtfs[species[i]].query('Feature=="exon"')
        gtf_exon_combined_temp.loc[:,'Chromosome'] = species[i] + '_' + gtf_exon_combined_temp.Chromosome.apply(lambda s:str(s))
        gtf_exon_combined = pd.concat([gtf_exon_combined,gtf_exon_combined_temp])
    gene_biotypes = gtf_exon_combined.Attributes.apply(lambda s: get_attribute(s,'gene_biotype'))
    #gtf_exon_combined = gtf_exon_combined.iloc[np.where(gene_biotypes.isin(gene_biotypes_to_keep).values)]
    gtf_exon_combined.index = range(len(gtf_exon_combined))
    gtf_exon_combined.to_csv(output_dir + '/exons.gtf',sep='\t',index=False)
    
    if not splicing:
        gtf_gene_combined = gtf_exon_combined.copy(deep=True)
        gtf_gene_combined['Feature'] = 'gene'
        gtf_gene_combined.to_csv(output_dir + '/genes.gtf',sep='\t',index=False)
    # Get locations of genes. We are using the longest possible span of different transcripts here
    gtf_gene_combined.loc[:,'gene_id'] = gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_id'))
    gene_starts = gtf_gene_combined.groupby('gene_id').Start.apply(min)
    gene_ends = gtf_gene_combined.groupby('gene_id').End.apply(max)
    chroms = gtf_gene_combined.groupby('gene_id').Chromosome.apply(lambda s:list(s)[0])
    strands = gtf_gene_combined.groupby('gene_id').Strand.apply(lambda s:list(s)[0])
    
    gtf_dict_stepsize = 10000
    # Create a dictionary for each "bin" of the genome, that maps to a list of genes within or overlapping
    # that bin. The bin size is determined by gtf_dict_stepsize.
    starts_rounded = gene_starts.apply(lambda s:np.floor(s/gtf_dict_stepsize)*gtf_dict_stepsize).values
    ends_rounded = gene_ends.apply(lambda s:np.ceil(s/gtf_dict_stepsize)*gtf_dict_stepsize).values
    gene_ids = gene_starts.index
    start_dict = gene_starts.to_dict()
    end_dict = gene_ends.to_dict()
    gene_dict = defaultdict(list)
    for i in range(len(gene_starts)):
        cur_chrom = chroms[i]
        cur_strand = strands[i]
        cur_start = int(starts_rounded[i])
        cur_end = int(ends_rounded[i])
        cur_gene_id = gene_ids[i]
        for coord in range(cur_start,cur_end+1,gtf_dict_stepsize):
            if not (cur_gene_id in gene_dict[cur_chrom + ':' +  str(coord)]):
                gene_dict[cur_chrom + ':' +  str(coord)+':'+cur_strand].append(cur_gene_id)
                
    # Create a dictionary from genes to exons
    exon_gene_ids = gtf_exon_combined.Attributes.apply(lambda s: get_attribute(s,'gene_id')).values
    exon_starts = gtf_exon_combined.Start.values
    exon_ends = gtf_exon_combined.End.values
    
    exon_gene_start_end_dict = defaultdict(dict)
    for i in range(len(exon_gene_ids)):
        cur_gene_id = exon_gene_ids[i]
        cur_exon_start = exon_starts[i]
        cur_exon_ends = exon_ends[i]
        exon_gene_start_end_dict[cur_gene_id][cur_exon_start] = cur_exon_ends
        
    gene_id_to_gene_names = dict(zip(gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_id')),
                                     gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_name'))))
    gene_id_to_genome = dict(zip(gtf_gene_combined.Attributes.apply(lambda s: get_attribute(s,'gene_id')),
                                 gtf_gene_combined.Chromosome.apply(lambda s:s.split('_')[0])))
    gene_id_to_strand = dict(zip(gtf_gene_combined.Attributes.apply(lambda s:get_attribute(s,'gene_id')).values,
                                 gtf_gene_combined.Strand.values))
    gene_id_to_chrom = dict(zip(gtf_gene_combined.Attributes.apply(lambda s:get_attribute(s,'gene_id')).values,
                                 gtf_gene_combined.Chromosome.values))
    gene_id_to_biotype = dict(zip(gtf_gene_combined.Attributes.apply(lambda s:get_attribute(s,'gene_id')).values,
                                  gtf_gene_combined.Attributes.apply(lambda s:get_attribute(s,'gene_biotype')).values))
    
    #Save dictionary with gene info
    gene_info = {'gene_bins':gene_dict,
                 'genes_to_exons':exon_gene_start_end_dict,
                 'gene_starts': start_dict,
                 'gene_ends': end_dict,
                 'gene_id_to_name': gene_id_to_gene_names,
                 'gene_id_to_genome':gene_id_to_genome,
                 'gene_id_to_chrom':gene_id_to_chrom,
                 'gene_id_to_strand':gene_id_to_strand,
                 'gene_id_to_biotype':gene_id_to_biotype
                }
    
    with open(output_dir+ '/gene_info.pkl', 'wb') as f:
        pickle.dump(gene_info, f, pickle.HIGHEST_PROTOCOL)
        
def generate_STAR_index(output_dir, nthreads,genomeSAindexNbases,splicing):
    splicing = (splicing=='True')
    if splicing:
        star_command = """STAR  --runMode genomeGenerate --genomeDir {0} --genomeFastaFiles {0}/genome.fa --sjdbGTFfile {0}/exons.gtf --runThreadN {1} --limitGenomeGenerateRAM 24000000000 --genomeSAindexNbases {2}""".format(output_dir, nthreads, genomeSAindexNbases)
    else:
        star_command = """STAR  --runMode genomeGenerate --genomeDir {0} --genomeFastaFiles {0}/genome.fa --runThreadN {1} --limitGenomeGenerateRAM 24000000000 --genomeSAindexNbases {2}""".format(output_dir, nthreads, genomeSAindexNbases)
    rc = subprocess.call(star_command, shell=True)
    return rc

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

def get_min_edit_dists(bc,edit_dict,max_d=3):
    """Returns a list of nearest edit dist seqs
    Input 8nt barcode, edit_dist_dictionary
    Output <list of nearest edit distance seqs>, <edit dist>"""
    bc_matches = edit_dict[0][bc]
    edit_dist = 0
    if (len(bc_matches)==0) and (max_d>=1):
        edit_dist+=1
        bc_matches = edit_dict[1][bc]
    if (len(bc_matches)==0) and (max_d>=2):
        edit_dist+=1
        bc_matches = edit_dict[2][bc]
    if (len(bc_matches)==0) and (max_d>=3):
        edit_dist+=1
        bc_matches = edit_dict[3][bc]
    return bc_matches,edit_dist

def preprocess_fastq(fastq1, fastq2, output_dir, chemistry='v1', bc_edit_dist=3, **params):
    """
    Performs all the steps before running the alignment. Temporary files
    saved in output_dir.
    """
    with open(PATH + '/barcodes/bc_dict_v1.pkl', 'rb') as f:
        edit_dict_v1 = pickle.load(f)
    with open(PATH + '/barcodes/bc_dict_v2.pkl', 'rb') as f:
        edit_dict_v2 = pickle.load(f)
    with open(PATH + '/barcodes/bc_dict_v3.pkl', 'rb') as f:
        edit_dict_v3 = pickle.load(f)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bc_edit_dist = int(bc_edit_dist)

    # Read in barcode sequences
    bc_8nt_v1 = pd.read_csv(PATH + '/barcodes/bc_8nt_v1.csv',names=['barcode'],index_col=0).barcode.values
    bc_8nt_v2 = pd.read_csv(PATH + '/barcodes/bc_8nt_v2.csv',names=['barcode'],index_col=0).barcode.values

    if chemistry=='v1':
        bc1_edit_dict = edit_dict_v1
        bc2_edit_dict = edit_dict_v1
        bc3_edit_dict = edit_dict_v1
        # Amplicon sequence
        amp_seq = 'NNNNNNNNNNIIIIIIIIGTGGCCGATGTTTCGCATCGGCGTACGACTIIIIIIIIATCCACGTGCTTGAGAGGCCAGAGCATTCGIIIIIIII'
    elif chemistry=='v2':
        bc1_edit_dict = edit_dict_v1
        bc2_edit_dict = edit_dict_v1
        bc3_edit_dict = edit_dict_v2

        # Amplicon sequence
        amp_seq = 'NNNNNNNNNNIIIIIIIIGTGGCCGATGTTTCGCATCGGCGTACGACTIIIIIIIIATCCACGTGCTTGAGACTGTGGIIIIIIII'

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

    print(bc_starts)

    def get_perfect_bc_counts(fastq2,n_reads=2000000,reads_in_cells_thresh=0.92):
        quality_scores = []
        seqs = []
        with gzip.open(fastq2) as f:
            for i in range(n_reads):
                f.readline()
                seq = f.readline().decode()[:-1]
                f.readline()
                qual = f.readline()
                seqs.append(seq)
                quality_scores.append(qual)
                if i %100000==0:
                    print(i,end=' ')
        seqs = pd.Series(seqs)
        bc_df = pd.DataFrame()
        bc_df['bc1'] = seqs.str.slice(bc_starts[0],bc_starts[0]+8)
        bc_df['bc2'] = seqs.str.slice(bc_starts[1],bc_starts[1]+8)
        bc_df['bc3'] = seqs.str.slice(bc_starts[2],bc_starts[2]+8)
        bc_df['bc1_valid'] = bc_df['bc1'].apply(lambda s: s in bc_8nt_v1)
        bc_df['bc2_valid'] = bc_df['bc2'].apply(lambda s: s in bc_8nt_v1)
        if chemistry=='v1':
            bc_df['bc3_valid'] = bc_df['bc3'].apply(lambda s: s in bc_8nt_v1)
        elif chemistry=='v2':
            bc_df['bc3_valid'] = bc_df['bc3'].apply(lambda s: s in bc_8nt_v2)

        counts = bc_df.query('bc1_valid & bc2_valid & bc3_valid')\
                      .groupby(['bc1','bc2','bc3']).size().sort_values(ascending=False)
        count_threshold = max(5, counts.iloc[abs(counts.cumsum()/counts.sum()-reads_in_cells_thresh).values.argmin()])
        return counts.to_dict(),count_threshold


    counts,count_threshold = get_perfect_bc_counts(fastq2,n_reads=2000000,reads_in_cells_thresh=0.92)

    def correct_barcodes(bc1,bc2,bc3, counts, count_thresh,
                     bc1_dict=bc1_edit_dict,
                     bc2_dict=bc2_edit_dict,
                     bc3_dict=bc3_edit_dict
                    ):
        bc1_matches,edit_dist1 = get_min_edit_dists(bc1,bc1_dict,max_d=bc_edit_dist)
        bc2_matches,edit_dist2  = get_min_edit_dists(bc2,bc2_dict,max_d=bc_edit_dist)
        bc3_matches,edit_dist3  = get_min_edit_dists(bc3,bc3_dict,max_d=bc_edit_dist)
        # Check if any barcode matches have a counts above the threshold

        if 0==edit_dist1==edit_dist2==edit_dist3:
            return bc1_matches[0],bc2_matches[0],bc3_matches[0]
        else:
            matches = 0
            for bc1_m in bc1_matches:
                for bc2_m in bc2_matches:
                    for bc3_m in bc3_matches:
                        try:
                            cur_counts = counts[(bc1_m,bc2_m,bc3_m)]
                        except:
                            cur_counts = 0
                        if cur_counts>count_thresh:
                            bc1_fixed = bc1_m
                            bc2_fixed = bc2_m
                            bc3_fixed = bc3_m
                            matches += 1
            #print('Fixed:',matches,\
            #      'Nearest_bcs',(len(bc1_matches),len(bc2_matches),len(bc3_matches)),\
            #      'Edit_dist',(edit_dist1,edit_dist2,edit_dist3))
            if matches==1:
                return (bc1_fixed,bc2_fixed,bc3_fixed)
            else:
                return '','',''

    fastq_reads = 0
    fastq_valid_barcode_reads = 0
    bc1_Q30_sum = 0
    bc2_Q30_sum = 0
    bc3_Q30_sum = 0
    umi_Q30_sum = 0
    cDNA_Q30_sum = 0
    with gzip.open(fastq1,'rb') as f1, gzip.open(fastq2,'rb') as f2, open(output_dir + '/single_cells_barcoded_head.fastq','w') as fout:
        while True:
            header2 = f2.readline()
            if len(header2)==0:
                break
            seq2 = f2.readline().decode("utf-8")
            bc1 = seq2[bc_starts[0]:bc_starts[0]+bc_len]
            bc2 = seq2[bc_starts[1]:bc_starts[1]+bc_len]
            bc3 = seq2[bc_starts[2]:bc_starts[2]+bc_len]
            umi = seq2[:10]
            strand2 = f2.readline()
            qual2 = f2.readline().decode("utf-8")
            bc1,bc2,bc3 = correct_barcodes(bc1,bc2,bc3, counts, count_threshold)
            cellbc_umi = bc3 + bc2 + bc1 +'_' + umi

            header1 = f1.readline().decode("utf-8")
            seq1 = f1.readline().decode("utf-8")
            strand1 = f1.readline().decode("utf-8")
            qual1 = f1.readline().decode("utf-8")
            
            if len(cellbc_umi)==35:
                TSO_location = seq1.find('AAGCAGTGGTATCAACGCAGAGTGAATGGG')
                if 0<=TSO_location<20:
                    seq1 = seq1[TSO_location+30:]
                    qual1 = qual1[TSO_location+30:]
                header1 = '@' + bc1 + bc2 + bc3 +'_' + umi + '_' + qual2[:10] + '_' + header1[1:]
                fout.write(header1)
                fout.write(seq1)
                fout.write(strand1)
                fout.write(qual1)
                fastq_valid_barcode_reads += 1
            fastq_reads += 1
            # bc1 refers to the first barcode seen in the sequencing read, but is actually
            # bc3 in terms of rounds
            if fastq_reads<1000000:
                bc1_Q30_sum += np.mean([ord(c)>62 for c in qual2[bc_starts[0]:bc_starts[0]+bc_len]])
                bc2_Q30_sum += np.mean([ord(c)>62 for c in qual2[bc_starts[1]:bc_starts[1]+bc_len]])
                bc3_Q30_sum += np.mean([ord(c)>62 for c in qual2[bc_starts[2]:bc_starts[2]+bc_len]])
                umi_Q30_sum += np.mean([ord(c)>62 for c in qual2[:10]])
                cDNA_Q30_sum += np.mean([ord(c)>62 for c in qual1[:-1]])
    with open(output_dir + '/sequencing_stats.txt', 'w') as f:
        f.write('bc1_Q30\t%0.4f\n' %(bc3_Q30_sum/min(fastq_reads,1000000)))
        f.write('bc2_Q30\t%0.4f\n' %(bc2_Q30_sum/min(fastq_reads,1000000)))
        f.write('bc3_Q30\t%0.4f\n' %(bc1_Q30_sum/min(fastq_reads,1000000)))
        f.write('umi_Q30\t%0.4f\n' %(umi_Q30_sum/min(fastq_reads,1000000)))
        f.write('cDNA_Q30\t%0.4f\n' %(cDNA_Q30_sum/min(fastq_reads,1000000)))
    with open(output_dir + '/pipeline_stats.txt', 'w') as f:
        f.write('fastq_reads\t%d\n' %fastq_reads)
        f.write('fastq_valid_barcode_reads\t%d\n' %fastq_valid_barcode_reads)

    return 0

def run_star(genome_dir, output_dir, nthreads):
    """ Align reads using STAR.
    """

    nthreads = int(nthreads)
    rc = subprocess.call(STAR_PATH + """ --genomeDir {0}/ --runThreadN {2} --readFilesIn {1}/single_cells_barcoded_head.fastq --outFileNamePrefix {1}/single_cells_barcoded_head""".format(genome_dir, output_dir, nthreads), shell=True)
    
    # Add alignment stats to pipeline_stats
    with open(output_dir + '/single_cells_barcoded_headLog.final.out') as f:
        for i in range(8):
            f.readline()
        unique_mapping = int(f.readline().split('\t')[1][:-1])
        for i in range(14):
            f.readline()
        multimapping = int(f.readline().split('\t')[1][:-1])
    with open(output_dir + '/pipeline_stats.txt', 'a') as f:
        f.write('uniquely_aligned\t%d\n' %unique_mapping)
        f.write('multimapping\t%d\n' %multimapping)

    return rc

def sort_sam(output_dir, nthreads):
    """ Sort samfile by header (cell_barcodes, umi) """
    nthreads = int(nthreads)
    rc = subprocess.call(SAMTOOLS_PATH + """ sort -n -@ {1} -T {0}/single_cells_barcoded_headAligned.sort -o {0}/single_cells_barcoded_headAligned.sorted.bam {0}/single_cells_barcoded_headAligned.out.sam""".format(output_dir, nthreads), shell=True)
    os.remove("""{0}/single_cells_barcoded_headAligned.out.sam""".format(output_dir))
    return rc
