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
              --nthreads 16 \
              --sample sample_name1 A1:B6 \
              --sample sample_name2 A7:B12 \
              --sample sample_name3 C1:D6 \
              --sample sample_name4 C7:D12
~~~~

## Merging Sublibraries into a Single Matrix

~~~~
split-seq combine --output_dir <output_dir> \
                  --sublibraries <path_to_sublibrary1> <path_to_sublibrary2> ...
                  --chemistry v2
                  --genome_dir <path_to_genome_dir>
                  --sample sample_name1 <wells>
~~~~

## Outputs

Running `split-seq all` with `--output_dir <output_dir>` generates three output folders: `<output_dir>`, `<output_dir>DGE_filtered`, and `<output_dir>DGE_unfiltered`

The first folder contains the read mappings and read assignments. Some important files:

- `read_assignments.csv` - this is a table that contains the gene assignment for every cell barcode-UMI combination.
- `single_cells_barcoded_head.fastq.gz` - this contains all reads in read1 that have a valid cell barcode. Each read is labeled with its barcode-UMI combination.
- `single_cells_barcoded_headAligned.sorted.bam` - this file contains the alignment to the genome for all reads in `single_cells_barcoded_head.fastq.gz`.

The `DGE_filtered` and `DGE_unfiltered` folders contain digital gene expression matrices. In `DGE_filtered`, the cells are filtered by a minimum read threshold, and only cells pasing that threshold are included.

In these two folders, `DGE.mtx` is a sparse matrix (Matrix Market format) of shape cells by genes that contains the gene expression of every gene for each cell. `genes.csv` contains the name of each gene.

## References
