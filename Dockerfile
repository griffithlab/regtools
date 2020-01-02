################################################################################
##################### Set Inital Image to work from ############################

# work from latest LTS ubuntu release
FROM ubuntu:18.04

# set variables
ENV r_version 3.6.0

# run update
RUN apt-get update -y && apt-get install -y \
  gfortran \
  libreadline-dev \
  libpcre3-dev \
  libcurl4-openssl-dev \
  build-essential \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  openjdk-8-jdk \
  wget \
  libssl-dev \
  libxml2-dev \
  libnss-sss \
  git \
  build-essential \
  cmake \
  python3

################################################################################
##################### Add Container Labels #####################################
LABEL "Regtools_License"="MIT"
LABEL "Description"="Software package which integrate DNA-seq and RNA-seq data\
                     to help interpret mutations in a regulatory and splicing\
                     context."

################################################################################
##################### Install Regtools #########################################

# clone git repository
RUN git clone https://github.com/griffithlab/regtools.git

# make a build directory for regtools
WORKDIR /regtools/
RUN mkdir build

# compile from source
RUN cd /regtools/build && cmake ..
RUN cd /regtools/build && make

################################################################################
###################### set environment path    #################################

# add regtools executable to path
ENV PATH="/regtools/build:${PATH}"

################################################################################
####################### Install R ##############################################

# change working dir
WORKDIR /usr/local/bin

# install R
RUN wget https://cran.r-project.org/src/base/R-3/R-${r_version}.tar.gz
RUN tar -zxvf R-${r_version}.tar.gz
WORKDIR /usr/local/bin/R-${r_version}
RUN ./configure --prefix=/usr/local/ --with-x=no
RUN make
RUN make install

# install R packages
RUN R --vanilla -e 'install.packages(c("data.table", "plyr", "tidyverse"), repos = "http://cran.us.r-project.org")'
