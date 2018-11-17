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

def preprocess_fastq(fastq1, fastq2, output_dir, picard_dir, **params):
    """
    Performs all the steps before running the alignment. Temporary files
    saved in output_dir.
    """
    picard_path = os.path.join(picard_dir, 'picard.jar')
    # calls picard to convert fastq to an unaligned bam file
    rc = subprocess.call('java -jar {3} FastqToSam F1={0} F2={1} O={2}/single_cells.bam SM=1'.format(fastq1, fastq2, output_dir, picard_path),
            shell=True)
    if rc == 1:
        return 'failed: step 1'
    # copy header - trim read 2
    rc = subprocess.call('samtools view -H {0}/single_cells.bam > {0}/single_cells.trimmed.sam'.format(output_dir),
            shell=True)
    if rc == 1:
        return 'failed: step 1'
    # Trim read2 to UMI_BC3_BC2_BC1 if bases 81-86=="CATTCG" (no dephasing):
    # TODO: should the awk script be modifiable?
    rc = subprocess.call("""
            samtools view {0}/single_cells.bam | awk 'NR%2==1{a=$0}NR%2==0 && substr($10,81,6)=="CATTCG"{print a"\n"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"substr($10,1,18)substr($10,49,8)substr($10,87,8)"\t"substr($11,1,18)substr($11,49,8)substr($11,87,8)"\t"$12}' >> {0}/single_cells.trimmed.sam
            """.format(output_dir), shell=True)
    # TODO: filter based on quality scores in barcode and UMI - allow only 1 base with phred<20 in each of BC and UMI.
    # Tag headers with BC3_BC2_BC1_UMI
    # Make reads unmated:
    # Convert sam to fastq
    return 0

def run_star(genome_dir, output_dir):
    """
    """
    rc = subprocess.call("""STAR --genomeDir {0}/ --runThreadN 8 --readFilesIn {1}/single_cells_barcoded_head.fastq --outFileNamePrefix {1}/single_cells_barcoded_head""".format(genome_dir, output_dir), shell=True)
    return rc

def run_postprocessing(input_dir, output_dir):
    """
    """
