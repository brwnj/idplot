FROM continuumio/miniconda3:4.8.2

LABEL version="1.1.2" maintainer="Joe Brown <brwnjm@gmail.com>"

RUN apt-get update && apt-get install -y g++ make procps
RUN conda create -n main -c bioconda -c conda-forge python==3.8.3 mafft==7.471 hyphy==2.5.15
RUN echo "source activate main" > ~/.bashrc
ENV PATH /opt/conda/envs/main/bin:$PATH

COPY ./docker-build.sh .
RUN bash docker-build.sh
