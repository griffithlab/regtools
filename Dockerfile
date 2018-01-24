################################################################################
##################### Set Inital Image to work from ############################
FROM ubuntu:16.04

################################################################################
##################### Add Container Labels #####################################
LABEL "Regtools_License"="MIT"
LABEL "Description"="Software package which integrate DNA-seq and RNA-seq data\
                     to help interpret mutations in a regulatory and splicing\
                     context."

################################################################################
##################### Install System Dependencies ##############################
RUN apt-get update
RUN apt-get install -y git
RUN apt-get install -y build-essential
RUN apt-get install -y cmake
RUN apt-get install -y zlib1g-dev

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
