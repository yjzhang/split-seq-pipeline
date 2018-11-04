# miscellaneous tools
import os
import shlex
import subprocess
import sys

PATH = './'


def fastq_to_bam(fastq1, fastq2, output_dir):
    """
    calls picard to convert the fastq file to an unaligned bam file.
    """
    picard_path = os.path.join(PATH, 'picard.jar')
    subprocess.call('java -jar {3} FastqToSam F1={0} F2={1} O={2}/single_cells.bam SM=1'.format(fastq1, fastq2, output_dir, picard_path))

def download_genome(genome_dir, ref='hg19'):
    """
    Downloads the hg19 reference genome...
    """

def run_star(genome_dir):
    """
    """
