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
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
~~~~

Download mouse reference gtf file:<br>
~~~~
wget ftp://ftp.ensembl.org/pub/release-93/gtf/mus_musculus/Mus_musculus.GRCm38.93.gtf
gunzip Homo_sapiens.GRCh38.93.gtf.gz
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
~~~~


## References
