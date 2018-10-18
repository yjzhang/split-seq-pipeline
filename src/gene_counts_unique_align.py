import pandas as pd
import numpy as np
import scipy
from pylab import *
import os
import pysam
from collections import Counter, OrderedDict
from itertools import izip, imap
import scipy.stats
import gzip
import json
import seaborn as sb
import sys
import operator
fsize=14


exon_tagged_file=sys.argv[1]
gene_tagged_file=sys.argv[2]

#Fix barcodes with 1mm
barcodes = pd.read_csv('/net/shendure/vol8/projects/BCBL/Mixing/ref/barcodes_96.csv',index_col=0,names=['barcodes']).barcodes.values
barcode_mm_dict = {}
bases = 'ACGT'
for bc in barcodes:
	for pos in range(8):
		for b in bases:
			cur_bc_mm = bc[:pos]+b+bc[pos+1:]
			barcode_mm_dict[cur_bc_mm] = bc
		
#barcode_mm_dict = dict(zip(barcodes,barcodes))

def fix_mismatches(seq):
	output_seq = ''
	for i in range(3):
		try:
			fixed_seq = barcode_mm_dict[seq[i*8:(i+1)*8]]
		except:
			fixed_seq = ''
		output_seq += fixed_seq
	if len(output_seq)==24:
		return output_seq+seq[24:]
	else:
		return ''

#Load reads
samfile_exon = pysam.Samfile(exon_tagged_file)
samfile_gene = pysam.Samfile(gene_tagged_file)
chrom_dict = pd.DataFrame(samfile_exon.header['SQ']).SN.to_dict()
read_mappings = []
read_barcodes = []

headers = []
i = 0
for read_exon in samfile_exon:
	read_gene = samfile_gene.next()
	if i%100000==0:
		print i,
	GE = ''	
	if (not read_exon.is_secondary) & (read_exon.mapping_quality==255):
		cur_read_map = chrom_dict[read_exon.tid]
		try:
			GE = dict(read_exon.tags)['GE']
		except:
			try:
				GE = dict(read_gene.tags)['GE']+'_INTRONIC'
			except:
				pass
		if len(GE)>0:
			read_mappings.append(GE+'_'+cur_read_map[:5])
			read_barcodes.append(read_exon.qname[:35])
	i+=1

samfile_exon.close()
samfile_gene.close()

read_barcodes = pd.Series(read_barcodes)
read_mappings = pd.Series(read_mappings)

# Fix barcodes with 1nt mismatch from a given barcode
read_barcodes_fixed = read_barcodes.apply(lambda s:fix_mismatches(s))
read_mappings = read_mappings[read_barcodes_fixed.apply(len)==35]
read_barcodes_fixed = read_barcodes_fixed[read_barcodes_fixed.apply(len)==35]

read_assignment = read_mappings.groupby(read_barcodes_fixed.str.slice(0,24)+read_barcodes_fixed.str.slice(24)).apply(lambda x:max(Counter(x).iteritems(), key=operator.itemgetter(1))[0])
#read_assignment = read_mappings.groupby(read_barcodes_fixed.str.slice(0,24)+read_barcodes_fixed.str.slice(24)).apply(map_single_species_reads)

#### Use starcode to collapse UMIs within 1 nt of each other
(read_barcodes_fixed.str.slice(0,24)+read_barcodes_fixed.str.slice(25)).to_csv('read_barcodes.txt',
																			   index=False,
																			   header=False)

os.system("cat read_barcodes.txt | ~/src/starcode/starcode -i /dev/stdin -d 1 -t 5 --non-redundant > read_barcodes_filtered.txt")

read_barcodes_filtered = pd.read_csv('read_barcodes_filtered.txt',names=['barcode_seq']).barcode_seq
read_barcodes_filtered = read_barcodes_filtered.str.slice(0,24)+'_'+read_barcodes_filtered.str.slice(24)
read_assignment_filtered = read_assignment.ix[read_barcodes_filtered.values]
read_assignment_filtered.to_csv('read_assignment.csv')

#read_assignment_filtered.to_csv('gene_read_counts.csv')
#read_counts = read_assignment_filtered.groupby(pd.Series(read_assignment_filtered.index).str.slice(0,24).values).apply(lambda x:Counter(x))#.unstack().fillna(0).to_sparse()
#read_counts.columns = ['multiple_alignment']+list(read_counts.columns[1:])
#read_counts.index.rename('cell_barcode',inplace=True)
#read_counts.to_csv('gene_read_counts.csv')
