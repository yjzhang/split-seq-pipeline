FROM python:3.7-slim
MAINTAINER Yue Zhang <yjzhang@cs.washington.edu>

WORKDIR /split_seq_pipeline

COPY install_dependencies.sh /split_seq_pipeline

# TODO: add files to path

RUN sh install_dependencies.sh /split-seq-reqs

COPY . /split_seq_pipeline

RUN pip install .

ENV PATH "/split-seq-reqs/bin:PATH"

ENTRYPOINT ["split-seq"]

CMD ["--help"]
