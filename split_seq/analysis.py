import pandas as pd
import scipy.io as sio
import numpy as np
import scipy.sparse
import scipy
import gzip
import collections
from collections import defaultdict
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

def generate_dge_matrix(df,read_cutoff=100):
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
    ax.legend(bbox_to_anchor=(0.08,1.0),fontsize=fsize-1,handletextpad=0.025)
    if fig is None:
        return 0
    else:
        return fig,ax

def plot_read_thresh(read_counts,fig=None,ax=None):
    window = 4
    sorted_read_counts = pd.Series(np.log10(read_counts.sort_values(ascending=False).values))
    x = np.log10(sorted_read_counts.groupby(sorted_read_counts).size()[::-1].cumsum())
    y = pd.Series(sorted_read_counts.groupby(sorted_read_counts).size()[::-1].index)
    threshold = int((pd.Series(pd.rolling_mean(y.diff().values/x.diff().values,window)).idxmin()-window/2.))
    read_threshold = read_counts.sort_values(ascending=False)[threshold]
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
        sub_wells = np.unique(sub_wells)
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
 
def generate_all_dge_reports(output_dir, genome_dir, chemistry, samples):
    if len(samples)>0:
        for i in range(len(samples)):
            sample_name = samples[i][0]
            sub_wells = parse_wells(samples[i][1])
            generate_single_dge_report(output_dir,genome_dir,chemistry,sample_name=sample_name,sub_wells=sub_wells)
    else:
        generate_single_dge_report(output_dir,genome_dir,chemistry)

def generate_single_dge_report(output_dir,genome_dir,chemistry,sample_name='',sub_wells=None):

    with open(genome_dir +'/gene_info.pkl', 'rb') as f:
        gene_info = pickle.load(f)

    gene_dict = gene_info['gene_bins']
    exon_gene_start_end_dict = gene_info['genes_to_exons']
    start_dict = gene_info['gene_starts']
    end_dict = gene_info['gene_ends']
    gene_id_to_name = gene_info['gene_id_to_name']
    gene_id_to_genome = gene_info['gene_id_to_genome']

    bc_8nt = pd.read_csv('/home/ubuntu/split-seq-pipeline/split_seq/barcodes/bc_8nt_%s.csv'  %chemistry,
                         index_col=0,
                         names=['barcode']).barcode
    bc_8nt_dict = dict(zip(bc_8nt.values,list(range(48))+list(range(48))))
    bc_8nt_randhex_dt_dict = dict(zip(bc_8nt.values,['dt']*48+['randhex']*48))

    df = pd.read_csv(output_dir + '/read_assignment.csv')
    total_reads = df.shape[0]
    df['rt_type'] = df.cell_barcode.apply(lambda s:bc_8nt_randhex_dt_dict[s[16:24]])
    df['cell_barcode'] = df.cell_barcode.apply(lambda s:s[:16]+'_'+str(bc_8nt_dict[s[16:24]]))
    df['well'] = df.cell_barcode.apply(lambda s: int(s.split('_')[-1]))
    
    if not (sub_wells is None):
        df = df.query('well in @sub_wells')

    read_counts = df.groupby('cell_barcode').size()

    read_counts = df.groupby('cell_barcode').size()
    fig,ax,read_thresh = plot_read_thresh(read_counts)

    digital_count_matrix,all_genes,barcodes = generate_dge_matrix(df,read_cutoff=100)

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
    cell_df['well'] = pd.Series(barcodes).apply(lambda s: s.split('_')[-1])
    
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

    stat_dict = {}
    with open(output_dir + '/pipeline_stats.txt') as f:
        while True:
            line = f.readline()[:-1]
            if len(line)==0:
                break
            k,v = line.split('\t')
            stat_dict[k] = int(v)
    stat_dict['Estimated Number of Cells'] = len(barcodes)
    stat_dict['Mean Reads per Cell'] = (stat_dict['fastq_reads'] * df.shape[0]/total_reads)/len(barcodes)
    stat_dict['Number of Reads'] = stat_dict['fastq_reads'] * df.shape[0]/total_reads
    stat_dict['Sequencing Saturation'] = 1-df.shape[0]/df.counts.sum()
    stat_dict['Valid Barcode Fraction'] = stat_dict['fastq_valid_barcode_reads']/stat_dict['fastq_reads']
    stat_dict['Reads Mapped to Transcriptome'] = stat_dict['mapped_to_transcriptome']/stat_dict['fastq_valid_barcode_reads']
    for s in species:
        stat_dict['%s Fraction Reads in Cells' %s] = digital_count_matrix[:,species_gene_inds[s]].sum()/\
                                                     df.query('genome=="%s"' %s, engine='python').shape[0]
        stat_dict['%s Median UMIs per Cell' %s] = np.median(species_umi_counts[s].iloc[np.where(species_assignments==s)])
        stat_dict['%s Median Genes per Cell' %s] = np.median(species_gene_counts[s].iloc[np.where(species_assignments==s)])
        stat_dict['%s Number of Cells Detected' %s] = sum(species_assignments==s)
        stat_dict['%s Exonic Fraction' %s] = df.loc[np.where(df.cell_barcode.isin(barcodes).values)].query('genome=="%s"' %s).exonic.mean()
        stat_dict['%s dT Fraction' %s] = (df.loc[np.where(df.cell_barcode.isin(barcodes).values)]\
                                           .query('genome=="%s"' %s)\
                                           .rt_type=='dt').mean()
    stat_dict['Fraction Reads in Cells'] = digital_count_matrix.sum()/df.shape[0]

    stat_catagories = ['Estimated Number of Cells']
    for s in species:
        stat_catagories.append('%s Number of Cells Detected' %s)
    for s in species:
        stat_catagories.append('%s Median UMIs per Cell' %s)
    for s in species:
        stat_catagories.append('%s Median Genes per Cell' %s)
    stat_catagories += ['Mean Reads per Cell',
                        'Number of Reads',
                        'Sequencing Saturation',
                        'Valid Barcode Fraction',
                        'Reads Mapped to Transcriptome',
                        'Fraction Reads in Cells'
                       ]
    for s in species:
        stat_catagories.append('%s Fraction Reads in Cells' %s)
                       
        
    for s in species:
        stat_catagories.append('%s Exonic Fraction' %s)
    for s in species:
        stat_catagories.append('%s dT Fraction' %s)

    # Subsample reads
    species_read_proportions = df.groupby('genome').size()/df.groupby('genome').size().sum()
    gene_counts_subsampled_df = {}
    umi_counts_subsampled_df = {}
    for s in df.genome.unique():
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
                                                        .loc[barcodes[np.where(species_assignments==s)]]).median()
            umi_counts_subsampled[subsample] = (species_df[sub_sampled_counts>0]
                                                        .groupby('cell_barcode')
                                                        .umi.size()
                                                        .loc[barcodes[np.where(species_assignments==s)]]).median()
        gene_counts_subsampled_df[s] = pd.Series(gene_counts_subsampled).fillna(0)
        umi_counts_subsampled_df[s] = pd.Series(umi_counts_subsampled).fillna(0)

    # Generate summary PDF
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0,0.5,0.45,0.5])
    h = 0.975
    c =0
    for k in stat_catagories:
        if c<9:
            text2write = k+' '*int((32-len(k)))+str(int(np.round(stat_dict[k])))
        else:
            text2write = k+' '*int((32-len(k)))+'%0.3f' %stat_dict[k]
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

    ax = fig.add_axes([0.0,0.05,0.6,0.4])
    for s in species:
        gene_counts_subsampled_df[s].plot(label=s,ax=ax)
    ax.legend()
    ax.set_title('Median Genes per Cell')
    ax.set_xlabel('Sequencing Reads per Cell')

    ax = fig.add_axes([0.75,0.05,0.6,0.4])
    for s in species:
        umi_counts_subsampled_df[s].plot(label=s,ax=ax)
    ax.legend()
    ax.set_title('Median UMIs per Cell')
    ax.set_xlabel('Sequencing Reads per Cell')
    fig.savefig(output_dir +'/' + sample_name + 'analysis_summary.pdf',bbox_inches='tight')
