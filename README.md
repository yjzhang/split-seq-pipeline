# Analysis tools for split-seq

## Requirements

Requires python 3.

Additional software needed:

- [STAR](https://github.com/alexdobin/STAR)
- [samtools](https://github.com/samtools/samtools)

To install all dependencies, try running `install_dependencies.sh`, which installs dependencies to `~/split_seq_reqs/`.

To install the package: run `pip install -e .` (might need sudo).

## Running the pipeline

To see all options, run `split-seq -h`.

To run the entire pipeline on a given dataset:

`split-seq all --fq1 <fastq-1> --fq2 <fastq-2> --output_dir <output-dir> --genome-dir <genome dir with STAR index> --gtf_file <path to gtf file>`

## Examples

### running `mkref`

The `mkref` mode is used to create an index with STAR. This requires 4 arguments: `--genome`, `--genes`, `--fasta`, and `--output_dir`.

`--genome` is a string or list of space-separated strings describing the species, ex. `mm10`, `hg19`, or `mm10 hg19`.

`--genes` is a gtf file, or a list of space-separated gtf files, one per species.

`--fasta` is a

`--output_dir` is where the output index will be written to.

Example command: `split-seq mkref --genome mm10 --genes Mus_musculus.GRCm38.71.gtf --fasta mm10_all_chroms.fa --output_dir mm10_genome`

### running the pipeline

Say that the sequencing produces two fastq files: `1.fastq` and `2.fastq`. 

`split-seq all --fq1 1.fastq --fq2 2.fastq --output_dir output_data --genome-dir mm10_genome --gtf_file <path to gtf file>`

## References
