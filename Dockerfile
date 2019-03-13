FROM python:3.7
MAINTAINER Yue Zhang <yjzhang@cs.washington.edu>

WORKDIR /split_seq_pipeline

COPY . /split_seq_pipeline

# TODO: add files to path

RUN sh install_dependencies.sh /split-seq-reqs

RUN pip install -e .

RUN cp split-seq /split-seq-reqs/bin

ENV PATH "/split-seq-reqs/bin:$PATH"

ENTRYPOINT ["split-seq"]

CMD ["--help"]

# example run command: docker run -v  /data/:/data/ ayuezhang27/split-seq-pipeline all --fq1 /data/split_seq_test/SRR6750041_1.fastq.gz --fq2 /data/split_seq_test/SRR6750041_2.fastq.gz --output_dir /data/split_seq_test/output_docker --nthreads 16 --chemistry v1 --genome_dir /data/reference_genomes/mm10/

