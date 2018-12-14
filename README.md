# Analysis tools for split-seq

## Requirements

Requires python 3.

Additional software needed:

- [STAR](https://github.com/alexdobin/STAR)
- [samtools](https://github.com/samtools/samtools)

To install all dependencies, try running `install_dependencies.sh`, which installs dependencies to `~/split_seq_reqs/`.

To install the package: run `pip install -e .` (might need sudo).

## Generating a reference genome
Download human reference genome
`wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
`gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`

Download human reference gtf file:
`wget ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz`
`gunzip Homo_sapiens.GRCh38.93.gtf.gz`

Generate split-seq reference:
split-seq mkref --genome hg38 mm10 \
                --fasta hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa \
                --output_dir /home/ubuntu/sc/ref/hg38_mm10/ 
                --genes hg38/Homo_sapiens.GRCh38.93.gtf mm10/Mus_musculus.GRCm38.93.gtf 
                --nthreads 16


## Running the pipeline

To see all options, run `split-seq -h`.

## References
