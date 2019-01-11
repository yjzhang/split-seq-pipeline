# Analysis tools for split-seq

## Requirements

Requires python 3.

Additional software needed:

- [STAR](https://github.com/alexdobin/STAR)
- [samtools](https://github.com/samtools/samtools)

To install all dependencies, try running `install_dependencies.sh`, which installs dependencies to `~/split_seq_reqs/`.

To install the package: run `pip install -e .` (might need sudo).

## Generating a reference genome
Download human reference genome<br>
~~~~
wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
~~~~
Download human reference gtf file:<br>
~~~~
wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz
gunzip Homo_sapiens.GRCh38.93.gtf.gz
~~~~
Download mouse reference genome<br>
~~~~
wget ftp://ftp.ensembl.org/pub/release-93/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
~~~~

Download mouse reference gtf file:<br>
~~~~
wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf.gz
gunzip Mus_musculus.GRCm38.93.gtf.gz
~~~~

Generate split-seq reference:
~~~~
split-seq mkref --genome hg38 mm10 \
                --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa Mus_musculus.GRCm38.dna.primary_assembly.fa \
                --genes Homo_sapiens.GRCh38.93.gtf Mus_musculus.GRCm38.93.gtf 
                --output_dir <ref_path>/hg38_mm10/ 
                --nthreads 16
~~~~


## Running the pipeline

To see all options, run `split-seq -h`.
~~~~
split-seq all --fq1 input_R1.fastq.gz \
              --fq2 input_R2.fastq.gz \
              --output_dir <output_dir> \
              --chemistry v2 \
              --genome_dir <path_to_ref>/hg38_mm10/ \
              --nthreads 16
              --sample sample_name1 A1:B6
              --sample sample_name2 A7:B12
              --sample sample_name3 C1:D6
              --sample sample_name4 C7:D12
~~~~


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
