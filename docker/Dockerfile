FROM continuumio/miniconda3:4.8.2

LABEL version="1.0.0" maintainer="Joe Brown <brwnjm@gmail.com>"

RUN apt-get update && apt-get install -y procps
RUN conda create -n env -c bioconda minimap2==2.17 python==3.8.3 numpy==1.18.5 click==7.1.2
RUN echo "source activate env" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
