FROM continuumio/miniconda3:4.8.2

LABEL version="1.1.3" maintainer="Joe Brown <brwnjm@gmail.com>"

RUN apt-get update && apt-get install -y g++ make procps
RUN conda create -n main -c bioconda -c conda-forge python==3.8.3 mafft==7.471 hyphy==2.5.15 fasttree==2.1.10
RUN echo "source activate main" > ~/.bashrc
ENV PATH /opt/conda/envs/main/bin:$PATH
