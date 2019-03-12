FROM python:3.7
MAINTAINER Yue Zhang <yjzhang@cs.washington.edu>

WORKDIR /split_seq_pipeline

COPY . /split_seq_pipeline

# TODO: add files to path

RUN sh install-dependencies.sh /split-seq-reqs

RUN pip install .

ENV PATH "/split-seq-reqs/bin:PATH"

ENTRYPOINT ["split-seq"]

CMD ["--help"]
