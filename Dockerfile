################################################################################
##################### Set Inital Image to work from ############################

# work from latest LTS ubuntu release
FROM ubuntu:20.04

# set variables
ENV r_version 3.6.0
ENV TZ=US/Chicago
ENV DEBIAN_FRONTEND noninteractive

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
  python3 \
  python3-pip
  
################################################################################
##################### Add Container Labels #####################################
LABEL "Regtools_License"="MIT"
LABEL "Description"="Software package which integrate DNA-seq and RNA-seq data\
                     to help interpret mutations in a regulatory and splicing\
                     context."

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
RUN R --vanilla -e 'install.packages(c("data.table", "plyr", "tidyverse", "optparse"), repos = "http://cran.us.r-project.org")'

################################################################################
##################### Install SpliceAI #########################################

RUN pip3 install spliceai
RUN pip3 install --upgrade tensorflow
RUN pip3 install keras==2.4.3

################################################################################
##################### Install other python libraries ###########################

RUN pip3 install dfply
RUN pip3 install pandas
RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install argparse

################################################################################
##################### Install Regtools #########################################


# removed this due to docker build pulling the correct branch already and the below command actually overwriting the desired branch to master
# clone git repository
ADD . /regtools

# change to regtools to build it 

WORKDIR /regtools

# compile from source
RUN ls
RUN mkdir build && cd build && cmake .. && make

################################################################################
################### Make scripts executable ####################################

WORKDIR /regtools/scripts

RUN chmod ugo+x *

################################################################################
###################### set environment path    #################################

# add regtools executable to path
ENV PATH="/regtools/build:/regtools/scripts:/usr/local/bin:/usr/local/bin/R-${r_version}:${PATH}"

