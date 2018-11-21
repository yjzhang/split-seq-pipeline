# miscellaneous tools
import os
import shlex
import subprocess
import sys

from tqdm import tqdm

import HTSeq
import pysam

PATH = './'

def trim_filter_read(input_fastq_read1, input_fastq_read2,
        trimmed_read2_output,
        tagged_read1_output, **kwargs):
    """
    trims the read to UMI_BC1_BC2_BC3, and filters the read by quality score.

    This should be done for read2.
    """
    fastq_1_file = HTSeq.FastqReader(input_fastq_read1)
    fastq_2_file = HTSeq.FastqReader(input_fastq_read2)
    output_file = open(trimmed_read2_output, 'wa')
    output_reads = []
    output_indices = []
    for index, reads in tqdm(enumerate(zip(fastq_1_file, fastq_2_file))):
        read0 = reads[0]
        read = reads[1]
        if read[80:86] == 'CATTCG':
            # check for UMI/barcode quality
            quality = True
            read_umi = read[:10]
            bc3 = read[10:18]
            bc2 = read[48:56]
            bc1 = read[86:94]
            if sum(read_umi.qual < 20) > 1:
                quality = False
            elif sum(bc3.qual < 20) + sum(bc2.qual < 20) + sum(bc1.qual < 20) > 1:
                quality = False
            if quality:
                # TODO: tag read1 with the barcode + UMI?
                new_read = HTSeq.SequenceWithQualities(read_umi.seq + bc3.seq + bc2.seq + bc1.seq,
                    read.name,
                    read_umi.qualstr + bc3.qualstr + bc2.qualstr + bc1.qualstr,
                    'phred')
                output_indices.append(index)
                output_reads.append(new_read)
                new_read.write_to_fastq(output_file)
    return output_indices


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
