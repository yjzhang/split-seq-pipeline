import pandas as pd
import scipy.io as sio
import scipy.interpolate
import numpy as np
import scipy.sparse
import scipy
import gzip
import subprocess
import collections
from collections import defaultdict, Counter
import scipy.sparse as sp_sparse
import warnings
import pickle
import os
#warnings.filterwarnings('ignore')
import matplotlib
matplotlib.use('Agg')
import pylab as plt
fsize=14

PATH = os.path.dirname(__file__)

rc_dict = dict(zip(list('NACGT'),list('NTGCA')))
def reverse_complement(seq):
    return ''.join([rc_dict[s] for s in seq][::-1])

def generate_dge_matrix(df,read_cutoff=10):
    reads_per_cell = df.groupby(df.cell_barcode).size()
    cells = reads_per_cell[reads_per_cell>3]
    all_genes = pd.Series(df.gene.unique()).sort_values()
    all_genes.index = range(len(all_genes))
    gene_dict = dict(zip(all_genes.values,range(len(all_genes.values))))
    cell_dict = dict(zip(cells.index.values,range(len(cells.index.values))))
    rows,cols,vals = [],[],[]
    for bc,g in zip(df.cell_barcode.values,df.gene.values):
        try:
            cell_dict[bc]
        except:
            pass
        else:
            rows.append(cell_dict[bc])
            cols.append(gene_dict[g])
            vals.append(1)
    rows.append(len(cell_dict)-1)
    cols.append(len(gene_dict)-1)
    vals.append(0)
    digital_count_matrix = scipy.sparse.csr_matrix((vals,(rows,cols)),dtype=np.float64)
    thresholded_cells = np.array(digital_count_matrix.sum(1)).flatten()>read_cutoff
    digital_count_matrix = digital_count_matrix[thresholded_cells,:]
    expressed_genes = np.array(digital_count_matrix.sum(0)).flatten()>0
    all_genes = pd.Series(all_genes[expressed_genes])
    digital_count_matrix = digital_count_matrix[:,expressed_genes]
    barcodes = cells.index.values[thresholded_cells]
    return digital_count_matrix,all_genes,barcodes

def barnyard(cell_data,tickstep=10000,s=4,lim=None,ax=None,fig=None):
    species = cell_data.columns[:2]
    colors = [(0.8941176470588236, 0.10196078431372549, 0.10980392156862745),
              (0.21568627450980393, 0.49411764705882355, 0.7215686274509804),
              'gray']
    #colors = list(sb.color_palette('Set1',n_colors=2)) + ['gray']
    #sb.set_style("white")
    #sb.set_style("ticks")
    
    if ax is None:
        fig = figure(figsize=(3,3))
        ax = fig.add_subplot(111)
    counts1 = cell_data.iloc[:,0]
    counts2 = cell_data.iloc[:,1]
    cell_type1 = counts1>(counts2*9)
    cell_type2 = counts2>(counts1*9)
    mixed_cells = ~(cell_type1|cell_type2)
    plt.scatter(counts1[mixed_cells],
            counts2[mixed_cells],
            color=colors[2],
            s=s,
            label=None)
    plt.scatter(counts1[cell_type2],
            counts2[cell_type2],
            color=colors[0],
            s=s,
            alpha=1,
            label=None)
    plt.scatter(counts1[cell_type1],
            counts2[cell_type1],
            color=colors[1],
            s=s,
            label=None)
    plt.scatter([],[],
            color=colors[0],
            s=10,
            label='%d %s (%0.1f'%(sum(cell_type2),species[1],100*float(sum(cell_type2))/len(cell_type2))+'%)',
           )
    plt.scatter([],[],
            color=colors[1],
            label='%d %s (%0.1f'%(sum(cell_type1),species[0],100*float(sum(cell_type1))/len(cell_type1))+'%)',
            s=10)
    plt.scatter([],[],
            color=colors[2],
            label='%d Mixed (%0.1f'%(sum(mixed_cells),100*float(sum(mixed_cells))/len(mixed_cells))+'%)',
            s=10)

    if lim==None:
        lim = int((counts1+counts2).max()*1.1)
    ax.set_xticks(plt.arange(0,lim,tickstep))
    ax.set_yticks(plt.arange(0,lim,tickstep))
    ax.set_xticklabels(plt.arange(0,lim,tickstep),rotation=90)
    ax.axis([-int(lim/30.),lim,-int(lim/30.),lim])
    ax.set_xlabel('%s UMI Counts' %species[0],fontsize=fsize)
    ax.set_ylabel('%s UMI Counts' %species[1],fontsize=fsize)
    ax.tick_params(labelsize=fsize)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.legend(fontsize=fsize-1,handletextpad=0.025)
    if fig is None:
        return 0
    else:
        return fig,ax

def get_read_threshold(read_counts):
    window = 20
    read_counts = read_counts[read_counts>3]
    x = np.log10(np.arange(1,len(read_counts)+1))
    y = np.log10(read_counts).values
    f = scipy.interpolate.interp1d(x, y,kind='linear')
    x_hat = np.linspace(x.min(),x.max(),500)
    y_hat = f(x_hat)
    y_hat = pd.Series(index=x_hat,data=y_hat)
    y_hat_prime = (-y_hat).diff(window).iloc[window:].values
    threshold = 10**y_hat.iloc[np.argmax(y_hat_prime)]*0.5
    return threshold

def plot_read_thresh(read_counts,fig=None,ax=None):
    window = 4
    read_threshold = get_read_threshold(read_counts[read_counts>2])
    threshold = len(read_counts[read_counts>read_threshold])
    median_umis = read_counts.sort_values(ascending=False)[:threshold].median()
    if ax is None:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
    ax.plot(range(len(read_counts)),
            (read_counts.sort_values(ascending=False)).values,
            color='lightgray',
            linewidth=2)
    ax.plot(range(threshold),
            (read_counts.sort_values(ascending=False)).values[:threshold],
            color='g',
            linewidth=0,marker='.')
    ax.set_xscale('log')
    ax.set_yscale('log')
    _ = ax.set_xlabel('# Barcodes (logscale)')
    _ = ax.set_ylabel('# UMIs (logscale)')
    ax.text(1,10,' n_cells: %d\n read cutoff: %d\n median_umis: %d' %(threshold,read_threshold,median_umis))
    if fig is None:
        return read_threshold
    else:
        return fig,ax,read_threshold
    
def parse_wells(s):
    wells = np.arange(48,dtype=int).reshape(4,12)
    try:
        blocks = s.split(',')
        row_letter_to_number = {'A':0,'B':1,'C':2,'D':3}
        sub_wells = []
        for b in blocks:
            if ':' in b:
                start,end = b.split(':')
                s_row = row_letter_to_number[start[:1]]
                s_col = int(start[1:])-1
                e_row = row_letter_to_number[end[:1]]
                e_col = int(end[1:])
                sub_wells += list(wells[s_row:e_row+1,s_col:e_col].flatten())
            elif '-' in b:
                start,end = b.split('-')
                s_row = row_letter_to_number[start[:1]]
                s_col = int(start[1:])-1
                e_row = row_letter_to_number[end[:1]]
                e_col = int(end[1:])-1
                sub_wells += list(np.arange(wells[s_row,s_col],wells[e_row,e_col]+1))
        sub_wells = list(np.unique(sub_wells))
    except:
        sub_wells = 'Failed'
    return sub_wells

def check_valid_samples(samples):
    if len(samples)>0:
        for i in range(len(samples)):
            sample_name = samples[i][0]
            sub_wells = parse_wells(samples[i][1])
            if str(sub_wells) == 'Failed':
                return False
    return True
 
def generate_all_dge_reports(output_dir, genome_dir, chemistry, samples, sublibraries=None):
    if len(samples)>0:
        for i in range(len(samples)):
            sample_name = samples[i][0]
            sub_wells = parse_wells(samples[i][1])
            generate_single_dge_report(output_dir,genome_dir,chemistry,sample_name=sample_name,sub_wells=sub_wells,sublibraries=sublibraries)
    else:
        generate_single_dge_report(output_dir,genome_dir,chemistry,sublibraries=sublibraries)
    
    # gzip fastq file to save space
    if (not ('single_cells_barcoded_head.fastq.gz' in os.listdir(output_dir))) and (sublibraries is None):
        gzip_command = """gzip {0}/single_cells_barcoded_head.fastq""".format(output_dir)
        rc = subprocess.call(gzip_command, shell=True)

def generate_single_dge_report(output_dir,genome_dir,chemistry,sample_name='',sub_wells=None, read_thresh=None, sublibraries=None):
    
    # Load gene_info dictionary to assign genes to reads
    with open(genome_dir +'/gene_info.pkl', 'rb') as f:
        gene_info = pickle.load(f)

    gene_dict = gene_info['gene_bins']
    exon_gene_start_end_dict = gene_info['genes_to_exons']
    start_dict = gene_info['gene_starts']
    end_dict = gene_info['gene_ends']
    gene_id_to_name = gene_info['gene_id_to_name']
    gene_id_to_genome = gene_info['gene_id_to_genome']

    # Load the RT barcodes to collapse randhex with dT primers
    bc_8nt = pd.read_csv(PATH + '/barcodes/bc_8nt_%s.csv'  %chemistry,
                         index_col=0,
                         names=['barcode']).barcode
    bc_to_well = dict(zip(bc_8nt.values,range(96)))
    bc_8nt_dict = dict(zip(bc_8nt.values,list(range(48))+list(range(48))))
    bc_8nt_randhex_dt_dict = dict(zip(bc_8nt.values,['dt']*48+['randhex']*48))
    
    # Load the Ligation barcodes:
    bc_8nt_lig = pd.read_csv(PATH + '/barcodes/bc_8nt_v1.csv',
                         index_col=0,
                         names=['barcode']).barcode
    bc_to_well_lig = dict(zip(bc_8nt_lig.values,range(96)))

    # Load the read_assignment file
    if sublibraries is None:
        df = pd.read_csv(output_dir + '/read_assignment.csv')
        total_reads = df.shape[0]
        df['rt_type'] = df.cell_barcode.apply(lambda s:bc_8nt_randhex_dt_dict[s[16:24]])
        df['cell_barcode'] = df.cell_barcode.apply(lambda s:s[:16]+'_'+str(bc_8nt_dict[s[16:24]]))
        df['well'] = df.cell_barcode.apply(lambda s: int(s.split('_')[-1]))
    else:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        df = []
        i = 1
        for cur_dir in sublibraries:
            df.append(pd.read_csv(cur_dir + '/read_assignment.csv'))
            df[-1]['rt_type'] = df[-1].cell_barcode.apply(lambda s:bc_8nt_randhex_dt_dict[s[16:24]])
            df[-1]['cell_barcode'] = df[-1].cell_barcode.apply(lambda s:s[:16]+'_'+str(bc_8nt_dict[s[16:24]]))
            df[-1]['well'] = df[-1].cell_barcode.apply(lambda s: int(s.split('_')[-1]))
            df[-1]['cell_barcode'] = df[-1]['cell_barcode'] + '_s' + str(i)
            i+=1
        df = pd.concat(df)
        total_reads = df.shape[0]
    
    
    # Check if performing analysis on a subset of wells
    if not (sub_wells is None):
        df = df.query('well in @sub_wells')
    else:
        sub_wells = list(range(48))
    read_counts = df.groupby('cell_barcode').size().sort_values(ascending=False)
    fig,ax,read_thresh = plot_read_thresh(read_counts)

    digital_count_matrix,all_genes,barcodes = generate_dge_matrix(df,read_cutoff=10)

    gene_df = pd.DataFrame()
    gene_df['gene_id'] = all_genes
    gene_df['gene_name'] = all_genes.apply(lambda s:gene_id_to_name[s])
    gene_df['genome'] = all_genes.apply(lambda s:gene_id_to_genome[s])

    species = df.genome.unique()
    species_genes = {}
    species_gene_inds = {}
    species_umi_counts = {}
    species_gene_counts = {}
    for s in species:
        species_genes[s] = all_genes[all_genes.apply(lambda s:gene_id_to_genome[s])==s]
        species_gene_inds[s] = np.where(all_genes.apply(lambda s:gene_id_to_genome[s])==s)[0]
        species_umi_counts[s] = pd.Series(index=barcodes,
                                      data=np.array(digital_count_matrix[:,species_gene_inds[s]].sum(1)).flatten())
        species_gene_counts[s] = pd.Series(index=barcodes,
                                      data=np.array((digital_count_matrix[:,species_gene_inds[s]]>0).sum(1)).flatten())
        
    species_umi_counts = pd.DataFrame(species_umi_counts)
    species_gene_counts = pd.DataFrame(species_gene_counts)
    species_assignments = pd.Series(['multiplet' for i in range(len(barcodes))])
    for s in species:
        species_assignments.loc[np.where((species_umi_counts[s]/species_umi_counts.sum(1))>0.9)] = s

    cell_df = pd.DataFrame()
    cell_df['cell_barcode'] = pd.Series(barcodes)
    cell_df['species'] = species_assignments.values
    cell_df['rnd1_well'] = pd.Series(barcodes).apply(lambda s: int(s.split('_')[1]))
    cell_df['rnd2_well'] = pd.Series(barcodes).apply(lambda s: bc_to_well_lig[s[8:16]])
    cell_df['rnd3_well'] = pd.Series(barcodes).apply(lambda s: bc_to_well_lig[s[:8]])
    cell_df['umi_count'] = np.array(digital_count_matrix.sum(1)).flatten()
    cell_df['umi_count_50dup'] = cell_df['umi_count'] * 0.5/(1-df.shape[0]/df.counts.sum())
    cell_df['gene_count'] = np.array((digital_count_matrix>0).sum(1)).flatten()
    
    if len(sample_name)>0:
        sample_name = sample_name +'_'
    
    # Write unfiltered matrix data
    if not os.path.exists(output_dir + sample_name + 'DGE_unfiltered/'):
        os.makedirs(output_dir + sample_name + 'DGE_unfiltered/')

    gene_df.to_csv(output_dir + sample_name + 'DGE_unfiltered/genes.csv')
    cell_df.to_csv(output_dir + sample_name + 'DGE_unfiltered/cell_metadata.csv',index=False)
    sio.mmwrite(output_dir + sample_name + 'DGE_unfiltered/DGE.mtx',digital_count_matrix)

    # Filter based on automatic cutoff
    valid_cells = np.where(np.array(digital_count_matrix.sum(1)).flatten()>read_thresh)[0]
    digital_count_matrix = digital_count_matrix[valid_cells]
    barcodes = barcodes[valid_cells]
    cell_df = cell_df.iloc[valid_cells]

    # Write filtered matrix data
    if not os.path.exists(output_dir + sample_name + 'DGE_filtered/'):
        os.makedirs(output_dir + sample_name + 'DGE_filtered/')
    gene_df.to_csv(output_dir + sample_name + 'DGE_filtered/genes.csv')
    cell_df.to_csv(output_dir + sample_name + 'DGE_filtered/cell_metadata.csv',index=False)
    sio.mmwrite(output_dir + sample_name + 'DGE_filtered/DGE.mtx',digital_count_matrix)
        
    digital_count_matrix,all_genes,barcodes = generate_dge_matrix(df,read_cutoff=read_thresh)

    species_genes = {}
    species_gene_inds = {}
    species_umi_counts = {}
    species_gene_counts = {}
    for s in species:
        species_genes[s] = all_genes[all_genes.apply(lambda s:gene_id_to_genome[s])==s]
        species_gene_inds[s] = np.where(all_genes.apply(lambda s:gene_id_to_genome[s])==s)[0]
        species_umi_counts[s] = pd.Series(index=barcodes,
                                      data=np.array(digital_count_matrix[:,species_gene_inds[s]].sum(1)).flatten())
        species_gene_counts[s] = pd.Series(index=barcodes,
                                      data=np.array((digital_count_matrix[:,species_gene_inds[s]]>0).sum(1)).flatten())
        
    species_umi_counts = pd.DataFrame(species_umi_counts)
    species_gene_counts = pd.DataFrame(species_gene_counts)
    species_assignments = pd.Series(['multiplet' for i in range(len(barcodes))])
    for s in species:
        species_assignments.loc[np.where((species_umi_counts[s]/species_umi_counts.sum(1))>0.9)] = s
    species = np.unique(species_assignments.values)
    species = species[species!='multiplet']
    
    if sublibraries is None:
        # Calculate rRNA Percentage:
        kmer_len = 30
        rrna_sense_kmer_dict = {}
        rrna_antisense_kmer_dict = {}
        with open(PATH + '/rRNA.fa') as f:
            while True:
                line = f.readline()[:-1]
                if len(line)==0:
                    break
                if line[0]!='>':
                    for i in range(len(line)-kmer_len):
                        kmer = line[i:i+kmer_len]
                        rrna_sense_kmer_dict[kmer] = 0
                    line = reverse_complement(line)
                    for i in range(len(line)-kmer_len):
                        kmer = line[i:i+kmer_len]
                        rrna_antisense_kmer_dict[kmer] = 0
        
        kmer_len = 30
        mt_rrna_sense_kmer_dict = {}
        mt_rrna_antisense_kmer_dict = {}
        with open(PATH + '/mt_rRNA.fa') as f:
            while True:
                line = f.readline()[:-1]
                if len(line)==0:
                    break
                if line[0]!='>':
                    for i in range(len(line)-kmer_len):
                        kmer = line[i:i+kmer_len]
                        mt_rrna_sense_kmer_dict[kmer] = 0
                    line = reverse_complement(line)
                    for i in range(len(line)-kmer_len):
                        kmer = line[i:i+kmer_len]
                        mt_rrna_antisense_kmer_dict[kmer] = 0

        def search_kmers(seq,kmer_dict):
            found = False
            for i in range(0,41,10):
                try:
                    kmer_dict[seq[i:i+kmer_len]]
                    found = True
                except:
                    pass
            return found
        
        fastqfile = output_dir + '/single_cells_barcoded_head.fastq'
        if os.path.isfile(fastqfile):
            well_counts = {'total_counts':Counter(),
                           'rRNA_sense_counts':Counter(),
                           'rRNA_antisense_counts':Counter(),
                           'mt_rRNA_sense_counts':Counter(),
                           'mt_rRNA_antisense_counts':Counter()}
            read_lengths = Counter()
            with open(fastqfile) as f:
                for i in range(1000000):
                    header = f.readline()
                    if len(header)==0:
                        break
                    seq = f.readline()[:-1]
                    f.readline()
                    f.readline()
                    well = bc_to_well[header[17:17+8]]
                    well_counts['total_counts'][well] += 1
                    read_lengths[len(seq)]+=1
                    if search_kmers(seq,rrna_sense_kmer_dict):
                        well_counts['rRNA_sense_counts'][well] += 1
                    if search_kmers(seq,rrna_antisense_kmer_dict):
                        well_counts['rRNA_antisense_counts'][well] += 1
                    if search_kmers(seq,mt_rrna_sense_kmer_dict):
                        well_counts['mt_rRNA_sense_counts'][well] += 1
                    if search_kmers(seq,mt_rrna_antisense_kmer_dict):
                        well_counts['mt_rRNA_antisense_counts'][well] += 1
        else:
            fastqfile = fastqfile + '.gz'
            well_counts = {'total_counts':Counter(),
                           'rRNA_sense_counts':Counter(),
                           'rRNA_antisense_counts':Counter(),
                           'mt_rRNA_sense_counts':Counter(),
                           'mt_rRNA_antisense_counts':Counter()}

            read_lengths = Counter()
            with gzip.open(fastqfile) as f:
                for i in range(1000000):
                    header = f.readline().decode()
                    if len(header)==0:
                        break
                    seq = f.readline().decode()[:-1]
                    f.readline()
                    f.readline()
                    well = bc_to_well[header[17:17+8]]
                    well_counts['total_counts'][well] += 1
                    read_lengths[len(seq)]+=1
                    if search_kmers(seq,rrna_sense_kmer_dict):
                        well_counts['rRNA_sense_counts'][well] += 1
                    if search_kmers(seq,rrna_antisense_kmer_dict):
                        well_counts['rRNA_antisense_counts'][well] += 1
                    if search_kmers(seq,mt_rrna_sense_kmer_dict):
                        well_counts['mt_rRNA_sense_counts'][well] += 1
                    if search_kmers(seq,mt_rrna_antisense_kmer_dict):
                        well_counts['mt_rRNA_antisense_counts'][well] += 1
        read_len = max(read_lengths.keys())
        read_len_trimmed = read_len - 30
        tso_fraction = read_lengths[read_len_trimmed]/sum(read_lengths.values())
        cols = ['rRNA_sense_counts','rRNA_antisense_counts','total_counts']
        well_rrna_counts = pd.DataFrame(well_counts)[cols].reindex(sub_wells+list(np.array(sub_wells)+48)).fillna(0)
        well_rrna_counts_dt = pd.DataFrame(well_rrna_counts).reindex(sub_wells)
        well_rrna_counts_randhex = pd.DataFrame(well_rrna_counts).reindex(list(np.array(sub_wells)+48))
        well_rrna_fraction = (well_rrna_counts.T/well_rrna_counts.total_counts).T.iloc[:,:2]
        well_rrna_fraction_dt = (well_rrna_counts_dt.T/well_rrna_counts_dt.total_counts).T.iloc[:,:2]
        well_rrna_fraction_randhex = (well_rrna_counts_randhex.T/well_rrna_counts_randhex.total_counts).T.iloc[:,:2]
        rrna_fraction = well_rrna_counts.sum(0).iloc[:2]/well_rrna_counts.sum(0).iloc[2]
        rrna_fraction_dt = well_rrna_counts_dt.sum(0).iloc[:2]/well_rrna_counts_dt.sum(0).iloc[2]
        rrna_fraction_randhex = well_rrna_counts_randhex.sum(0).iloc[:2]/well_rrna_counts_randhex.sum(0).iloc[2]
        
        cols = ['mt_rRNA_sense_counts','mt_rRNA_antisense_counts','total_counts']
        well_mt_rrna_counts = pd.DataFrame(well_counts).fillna(0)[cols].reindex(sub_wells+list(np.array(sub_wells)+48))
        well_mt_rrna_counts_dt = pd.DataFrame(well_mt_rrna_counts).reindex(sub_wells)
        well_mt_rrna_counts_randhex = pd.DataFrame(well_mt_rrna_counts).reindex(list(np.array(sub_wells)+48))
        well_mt_rrna_fraction = (well_mt_rrna_counts.T/well_mt_rrna_counts.total_counts).T.iloc[:,:2]
        well_mt_rrna_fraction_dt = (well_mt_rrna_counts_dt.T/well_mt_rrna_counts_dt.total_counts).T.iloc[:,:2]
        well_mt_rrna_fraction_randhex = (well_mt_rrna_counts_randhex.T/well_mt_rrna_counts_randhex.total_counts).T.iloc[:,:2]
        mt_rrna_fraction = well_mt_rrna_counts.sum(0).iloc[:2]/well_mt_rrna_counts.sum(0).iloc[2]
        mt_rrna_fraction_dt = well_mt_rrna_counts_dt.sum(0).iloc[:2]/well_mt_rrna_counts_dt.sum(0).iloc[2]
        mt_rrna_fraction_randhex = well_mt_rrna_counts_randhex.sum(0).iloc[:2]/well_mt_rrna_counts_randhex.sum(0).iloc[2]
        pd.DataFrame(well_counts).fillna(0).reindex(sub_wells+list(np.array(sub_wells)+48)).to_csv(output_dir + sample_name + 'rRNA_count.csv')
        
        stat_dict = {}
        with open(output_dir + '/pipeline_stats.txt') as f:
            while True:
                line = f.readline()[:-1]
                if len(line)==0:
                    break
                k,v = line.split('\t')
                stat_dict[k] = int(v)
        with open(output_dir + '/sequencing_stats.txt') as f:
            while True:
                line = f.readline()[:-1]
                if len(line)==0:
                    break
                k,v = line.split('\t')
                stat_dict[k] = float(v)
        stat_dict['Estimated Number of Cells'] = len(barcodes)
        stat_dict['Mean Reads/Cell'] = (stat_dict['fastq_reads'] * df.shape[0]/total_reads)/len(barcodes)
        stat_dict['Number of Reads'] = stat_dict['fastq_reads'] * df.shape[0]/total_reads
        stat_dict['Sequencing Saturation'] = 1-df.shape[0]/df.counts.sum()
        stat_dict['Valid Barcode Fraction'] = stat_dict['fastq_valid_barcode_reads']/stat_dict['fastq_reads']
        stat_dict['BC1 (RT) >Q30'] = stat_dict['bc1_Q30']
        stat_dict['BC2 >Q30'] = stat_dict['bc2_Q30']
        stat_dict['BC3 >Q30'] = stat_dict['bc3_Q30']
        stat_dict['UMI >Q30'] = stat_dict['umi_Q30']
        stat_dict['cDNA >Q30'] = stat_dict['cDNA_Q30']
        stat_dict['Reads Mapped to rRNA'] = rrna_fraction.iloc[:2].sum()
        stat_dict['Reads Mapped to rRNA (dT)'] = rrna_fraction_dt.iloc[:2].sum()
        stat_dict['Reads Mapped to rRNA (randhex)'] = rrna_fraction_randhex.iloc[:2].sum()
        stat_dict['Reads Mapped to mt-rRNA'] = mt_rrna_fraction.iloc[:2].sum()
        stat_dict['Reads Mapped to mt-rRNA (dT)'] = mt_rrna_fraction_dt.iloc[:2].sum()
        stat_dict['Reads Mapped to mt-rRNA (randhex)'] = mt_rrna_fraction_randhex.iloc[:2].sum()
        stat_dict['TSO Fraction in Read1'] = tso_fraction
        stat_dict['Reads Mapped to Transcriptome'] = stat_dict['mapped_to_transcriptome']/stat_dict['fastq_valid_barcode_reads']
        for s in species:
            stat_dict['%s Fraction Reads in Cells' %s] = digital_count_matrix[:,species_gene_inds[s]].sum()/\
                                                         df.query('genome=="%s"' %s, engine='python').shape[0]
            stat_dict['%s Median UMIs/Cell' %s] = np.median(species_umi_counts[s].iloc[np.where(species_assignments==s)])
            stat_dict['%s Median UMIs/Cell @50%% Dup' %s]  = stat_dict['%s Median UMIs/Cell' %s] * 0.5 /stat_dict['Sequencing Saturation']
            stat_dict['%s Median Genes/Cell' %s] = np.median(species_gene_counts[s].iloc[np.where(species_assignments==s)])
            stat_dict['%s Number of Cells Detected' %s] = sum(species_assignments==s)
            stat_dict['%s Exonic Fraction' %s] = df.loc[np.where(df.cell_barcode.isin(barcodes).values)].query('genome=="%s"' %s).exonic.mean()
            stat_dict['%s dT Fraction' %s] = (df.loc[np.where(df.cell_barcode.isin(barcodes).values)]\
                                               .query('genome=="%s"' %s)\
                                               .rt_type=='dt').mean()
        stat_dict['Fraction Reads in Cells'] = digital_count_matrix.sum()/df.shape[0]
    else:
        stat_df = pd.concat([pd.read_csv(s +'/'+ sample_name + 'analysis_summary.csv',names=[s.split('/')[-2]]).T for s in sublibraries])
        num_reads_col = np.where(stat_df.columns=='Number of Reads')[0][0]
        stat_dict = ((stat_df.iloc[:,num_reads_col+1:].T*stat_df.iloc[:,num_reads_col]).sum(1)/stat_df.iloc[:,num_reads_col].sum()).T.to_dict()
        stat_dict['Estimated Number of Cells'] = len(barcodes)
        stat_dict['Mean Reads/Cell'] = (stat_df['Number of Reads'].sum() * df.shape[0]/total_reads)/len(barcodes)
        stat_dict['Number of Reads'] = stat_df['Number of Reads'].sum()
        for s in species:
            stat_dict['%s Fraction Reads in Cells' %s] = digital_count_matrix[:,species_gene_inds[s]].sum()/\
                                                                             df.query('genome=="%s"' %s, engine='python').shape[0]
            stat_dict['%s Median UMIs/Cell' %s] = np.median(species_umi_counts[s].iloc[np.where(species_assignments==s)])
            stat_dict['%s Median UMIs/Cell @50%% Dup' %s]  = stat_dict['%s Median UMIs/Cell' %s] * 0.5 /stat_dict['Sequencing Saturation']
            stat_dict['%s Median Genes/Cell' %s] = np.median(species_gene_counts[s].iloc[np.where(species_assignments==s)])
            stat_dict['%s Number of Cells Detected' %s] = sum(species_assignments==s)
        stat_dict['Fraction Reads in Cells'] = digital_count_matrix.sum()/df.shape[0]

    stat_catagories = ['Estimated Number of Cells']
    for s in species:
        stat_catagories.append('%s Number of Cells Detected' %s)
    for s in species:
        stat_catagories.append('%s Median UMIs/Cell' %s)    
    for s in species:
        stat_catagories.append('%s Median UMIs/Cell @50%% Dup' %s)
    for s in species:
        stat_catagories.append('%s Median Genes/Cell' %s)
    stat_catagories += ['Mean Reads/Cell',
                        'Number of Reads',
                        'Sequencing Saturation',
                        'Valid Barcode Fraction',
                        'BC1 (RT) >Q30',
                        'BC2 >Q30',
                        'BC3 >Q30',
                        'UMI >Q30',
                        'cDNA >Q30',
                        'Reads Mapped to rRNA',
                        'Reads Mapped to rRNA (dT)',
                        'Reads Mapped to rRNA (randhex)',
                        'Reads Mapped to mt-rRNA',
                        'Reads Mapped to mt-rRNA (dT)',
                        'Reads Mapped to mt-rRNA (randhex)',
                        'TSO Fraction in Read1',
                        'Reads Mapped to Transcriptome',
                        'Fraction Reads in Cells'
                       ]
    for s in species:
        stat_catagories.append('%s Fraction Reads in Cells' %s)
                       
        
    for s in species:
        stat_catagories.append('%s Exonic Fraction' %s)
    for s in species:
        stat_catagories.append('%s dT Fraction' %s)
    
    # Save summary stats to csv
    pd.Series(stat_dict).loc[stat_catagories].to_csv(output_dir + '/' + sample_name +  'analysis_summary.csv')
    
    # Subsample reads
    species_read_proportions = df.groupby('genome').size()/df.groupby('genome').size().sum()
    gene_counts_subsampled_df = {}
    umi_counts_subsampled_df = {}
    for s in species:
        seq_depth = species_read_proportions[s] * \
                    stat_dict['Number of Reads']/stat_dict['%s Number of Cells Detected' %s]
        gene_counts_subsampled = {0:0}
        umi_counts_subsampled = {0:0}
        subsample_depths = np.array(list(range(0,int(seq_depth),10000)) + [seq_depth],dtype=int)
        species_df = df.query('genome=="%s"'%s)
        for i in range(1, len(subsample_depths)):
            subsample = subsample_depths[i]
            subsample_fraction = subsample/seq_depth
            sub_sampled_counts = np.random.binomial(species_df.counts.values,subsample_fraction)
            gene_counts_subsampled[subsample] = (species_df[sub_sampled_counts>0]
                                                        .groupby('cell_barcode')
                                                        .gene.apply(lambda x:len(np.unique(x)))
                                                        .reindex(barcodes[np.where(species_assignments==s)])).median()
            umi_counts_subsampled[subsample] = (species_df[sub_sampled_counts>0]
                                                        .groupby('cell_barcode')
                                                        .umi.size()
                                                        .reindex(barcodes[np.where(species_assignments==s)])).median()
        gene_counts_subsampled_df[s] = pd.Series(gene_counts_subsampled).fillna(0)
        umi_counts_subsampled_df[s] = pd.Series(umi_counts_subsampled).fillna(0)
        
    # Get UMIs by well:
    umi_counts_by_rnd1_well = cell_df.groupby('rnd1_well').umi_count.median().reindex(range(48)).fillna(0)
    umi_counts_by_rnd2_well = cell_df.groupby('rnd2_well').umi_count.median().reindex(range(96)).fillna(0)
    umi_counts_by_rnd3_well = cell_df.groupby('rnd3_well').umi_count.median().reindex(range(96)).fillna(0)
    cell_counts_by_rnd1_well = cell_df.groupby('rnd1_well').size().reindex(range(48)).fillna(0)
    cell_counts_by_rnd2_well = cell_df.groupby('rnd2_well').size().reindex(range(96)).fillna(0)
    cell_counts_by_rnd3_well = cell_df.groupby('rnd3_well').size().reindex(range(96)).fillna(0)

    # Generate summary PDF
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0,0.5,0.45,0.5])
    h = 1
    c =0
    for k in stat_catagories:
        if c < (4*len(species)+3):
            text2write = k+' '*int((34-len(k)))+str(int(np.round(stat_dict[k])))
        else:
            text2write = k+' '*int((34-len(k)))+'%0.3f' %stat_dict[k]
        ax.text(-0.13,h,text2write,fontdict={'family':'monospace'},fontsize=11)
        h-=0.055
        c+=1
    ax.set_axis_off()
    ax = fig.add_axes([0.5,0.65,0.35,0.35])
    _ = plot_read_thresh(read_counts,ax=ax)
    ax.set_title(sample_name[:-1])

    if len(species)==2:
        ax = fig.add_axes([1,0.65,0.35,0.35])
        _ = barnyard(species_umi_counts,ax=ax)

    ax = fig.add_axes([0.5,0.1,0.35,0.35])
    for s in species:
        gene_counts_subsampled_df[s].plot(label=s,ax=ax)
    ax.legend()
    ax.set_title('Median Genes per Cell')
    ax.set_xlabel('Sequencing Reads per Cell')

    ax = fig.add_axes([1,0.1,0.35,0.35])
    for s in species:
        umi_counts_subsampled_df[s].plot(label=s,ax=ax)
    ax.legend()
    ax.set_title('Median UMIs per Cell')
    ax.set_xlabel('Sequencing Reads per Cell')
    
    fig.add_axes([0,-0.16,0.35,0.2])
    cm = plt.imshow(umi_counts_by_rnd1_well.values.reshape(4,12),cmap=plt.cm.Reds)
    plt.gca().set_xticks(range(0,12))
    plt.gca().set_xticklabels(range(1,13))
    plt.gca().set_yticks(range(0,4))
    plt.gca().set_yticklabels(list('ABCD'))
    plt.gca().set_title('Round1: Median UMIs per Well')
    plt.colorbar(cm,pad=0.02,aspect=10,shrink=0.7)
    
    fig.add_axes([0.5,-0.2,0.35,0.2])
    cm = plt.imshow(umi_counts_by_rnd2_well.values.reshape(8,12),cmap=plt.cm.Reds)
    plt.gca().set_xticks(range(0,12))
    plt.gca().set_xticklabels(range(1,13))
    plt.gca().set_yticks(range(0,8))
    plt.gca().set_yticklabels(list('ABCDEFGH'))
    plt.gca().set_title('Round2: Median UMIs per Well')
    plt.colorbar(cm,pad=0.02,aspect=10,shrink=0.7)
    
    fig.add_axes([1,-0.2,0.35,0.2])
    cm = plt.imshow(umi_counts_by_rnd3_well.values.reshape(8,12),cmap=plt.cm.Reds)
    plt.gca().set_xticks(range(0,12))
    plt.gca().set_xticklabels(range(1,13))
    plt.gca().set_yticks(range(0,8))
    plt.gca().set_yticklabels(list('ABCDEFGH'))
    plt.gca().set_title('Round3: Median UMIs per Well')
    plt.colorbar(cm,pad=0.02,aspect=10,shrink=0.7)
    
    fig.add_axes([0,-0.43,0.35,0.2])
    cm = plt.imshow(cell_counts_by_rnd1_well.values.reshape(4,12),cmap=plt.cm.Reds)
    plt.gca().set_xticks(range(0,12))
    plt.gca().set_xticklabels(range(1,13))
    plt.gca().set_yticks(range(0,4))
    plt.gca().set_yticklabels(list('ABCD'))
    plt.gca().set_title('Round1: Cells per Well')
    plt.colorbar(cm,pad=0.02,aspect=10,shrink=0.7)
    
    fig.add_axes([0.5,-0.47,0.35,0.2])
    cm = plt.imshow(cell_counts_by_rnd2_well.values.reshape(8,12),cmap=plt.cm.Reds)
    plt.gca().set_xticks(range(0,12))
    plt.gca().set_xticklabels(range(1,13))
    plt.gca().set_yticks(range(0,8))
    plt.gca().set_yticklabels(list('ABCDEFGH'))
    plt.gca().set_title('Round2: Cells per Well')
    plt.colorbar(cm,pad=0.02,aspect=10,shrink=0.7)
    
    fig.add_axes([1,-0.47,0.35,0.2])
    cm = plt.imshow(cell_counts_by_rnd3_well.values.reshape(8,12),cmap=plt.cm.Reds)
    plt.gca().set_xticks(range(0,12))
    plt.gca().set_xticklabels(range(1,13))
    plt.gca().set_yticks(range(0,8))
    plt.gca().set_yticklabels(list('ABCDEFGH'))
    plt.gca().set_title('Round3: Cells per Well')
    plt.colorbar(cm,pad=0.02,aspect=10,shrink=0.7)
    
    fig.savefig(output_dir +'/' + sample_name + 'analysis_summary.pdf',bbox_inches='tight')
