FROM ubuntu
LABEL maintainer="blodhi@korea.ac.kr" description="learning the docker environment"
USER root
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
apt-get install -y wget \
vim \
python \
htop \
python-pip \
python-dev \
unzip \
software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' 
RUN apt-get update && \
apt-get install -y r-base && \
apt-get -y upgrade
RUN apt-get install -y libcurl4-openssl-dev \
libssl-dev \
libxml2-dev \
libudunits2-dev \
git
RUN Rscript -e  'install.packages("BiocManager")' && \
Rscript -e 'BiocManager::install("ChIPseeker")'
RUN Rscript -e 'BiocManager::install("ChIPpeakAnno")' && \
Rscript -e 'BiocManager::install("SPIA")' && \
Rscript -e 'BiocManager::install("ReactomePA")'
Run Rscript -e 'BiocManager::install("optparse")' && \
Rscript -e 'BiocManager::install("clusterProfiler")'
Run Rscript -e 'BiocManager::install("org.Ce.eg.db")' && \
Rscript -e 'BiocManager::install("BSgenome.Celegans.UCSC.ce10")' && \
Rscript -e 'BiocManager::install("GenomicFeatures")'
RUN apt-get install -y libmariadb-client-lgpl
-dev  && \
Rscript -e 'BiocManager::install("RMariaDB")'