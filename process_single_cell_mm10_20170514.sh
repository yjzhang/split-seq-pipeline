#!/usr/bin/env bash

scripts_dir='/net/shendure/vol8/projects/BCBL/Mixing/bin'
dropseq_scripts_dir='/net/shendure/vol1/home/abros/bin/Drop-seq_tools-1.12'
dropseq_ref_dir='/net/shendure/vol8/projects/BCBL/Mixing/ref/mus10_ref/mm10'

fastq1=$1
fastq2=$2
output_dir=$3

mkdir $output_dir

# Convert fastq to unaligned bam:
java -jar ${dropseq_scripts_dir}/3rdParty/picard/picard.jar FastqToSam F1=$fastq1 F2=$fastq2 O=${output_dir}/single_cells.bam SM=1

cd ${output_dir}

# Copy header
echo "Trimming Read2"
samtools view -H single_cells.bam > single_cells.trimmed.sam

# Trim read2 to UMI_BC3_BC2_BC1 if bases 81-86=="CATTCG" (no dephasing):
samtools view single_cells.bam | awk 'NR%2==1{a=$0}NR%2==0 && substr($10,81,6)=="CATTCG"{print a"\n"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"substr($10,1,18)substr($10,49,8)substr($10,87,8)"\t"substr($11,1,18)substr($11,49,8)substr($11,87,8)"\t"$12}' >> single_cells.trimmed.sam

# Filter based on quality scores in barcode and UMI (allow only 1 base with phred<20 in each BC and UMI):
echo "Filtering UMIs and BCs based on quality scores"
python ${scripts_dir}/filter_qual_scores_umi.py single_cells.trimmed.sam single_cells_qs_filt.sam
rm single_cells.trimmed.sam

# Tag headers with BC3_BC2_BC1_UMI
echo "Tagging headers with BCs and UMIs"
samtools view -HS single_cells_qs_filt.sam > single_cells_barcoded_head.sam
samtools view -S single_cells_qs_filt.sam | awk 'NR%2==1{a=$0}NR%2==0 {print substr($10,11,24)"_"substr($10,1,10)"_"a}' >> single_cells_barcoded_head.sam
rm single_cells_qs_filt.sam

# Make reads unmated:
echo "Unmating reads"
sed -i -e 's/\t77/\t4/g' single_cells_barcoded_head.sam

# Convert sam to fastq
java -jar ${dropseq_scripts_dir}/3rdParty/picard/picard.jar SamToFastq INPUT=single_cells_barcoded_head.sam FASTQ=single_cells_barcoded_head.fastq
rm single_cells_barcoded_head.sam

# Run the STAR alignment
STAR --genomeDir ${dropseq_ref_dir}/ --runThreadN 8 --readFilesIn single_cells_barcoded_head.fastq --outFileNamePrefix single_cells_barcoded_head
rm single_cells_barcoded_head.fastq

# Sort the aligned bam file
java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar ${dropseq_scripts_dir}/3rdParty/picard/picard.jar SortSam INPUT=single_cells_barcoded_headAligned.out.sam OUTPUT=single_cells_barcoded_headAligned.sorted.bam SORT_ORDER=queryname
rm single_cells_barcoded_headAligned.out.sam

# Tag the reads with the gene and read type (only exonic get tagged with gene)
${dropseq_scripts_dir}/TagReadWithGeneExon O=star_gene_tagged.bam ANNOTATIONS_FILE=${dropseq_ref_dir}/mm10.refFlat TAG=GE CREATE_INDEX=true INPUT=single_cells_barcoded_headAligned.sorted.bam

# Tag the reads with the gene and read type (intronic also get tagged with gene)
${dropseq_scripts_dir}/TagReadWithGeneExon O=star_gene_tagged_full_gene.bam ANNOTATIONS_FILE=${dropseq_ref_dir}/mm10_full_genes.refFlat TAG=GE CREATE_INDEX=true INPUT=single_cells_barcoded_headAligned.sorted.bam
rm single_cells_barcoded_headAligned.sorted.bam

# Get the gene counts
python ${scripts_dir}/gene_counts_unique_align.py star_gene_tagged.bam star_gene_tagged_full_gene.bam

# Make the digital count matrix
python ${scripts_dir}/make_digital_count_matrix.py digital_gene_counts.mat ./

# Count the different read types (exon,intron,utr,intergenic)
python ${scripts_dir}/count_read_types.py star_gene_tagged.bam
