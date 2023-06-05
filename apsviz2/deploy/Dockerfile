##############
# Docker file for the creation of the adcirc support tools (AST).
##############
FROM continuumio/miniconda3

# get some credit
LABEL maintainer="jtilson@renci.org"

# update sources list, and install basic apps, one per line for better caching
RUN apt-get clean && apt-get --allow-releaseinfo-change update &&\
    apt-get install -qy vim

# update conda and set the download channel
RUN conda update conda && \
    conda config --add channels conda-forge

# tell conda what the shell is
RUN conda init bash

## create the non-root user
RUN useradd -m -d /home/nru -u 1000 nru

## make sure everything is read/write
RUN chmod 777 -R /home/nru

# make a directory for the repos and go there
RUN mkdir /home/nru/repo
WORKDIR /home/nru/repo

# get the AST repo - needs to be made public first
RUN git clone https://github.com/RENCI/AST.git
RUN git clone https://github.com/RENCI/EDSASTAPPS.git
RUN git clone https://github.com/RENCI/buoypy.git

# move to the deploy directory
WORKDIR /home/nru/repo/EDSASTAPPS/apsviz2/deploy

# Create the environment, initialize conda in bash config files, and install repos:
COPY environment.yml .
RUN conda env create -f environment.yml
                                
# install the pip requirements
RUN pip install -r pip_requirements.txt

# make sure everything is read/write
RUN chmod 777 -R /home/nru
                        
# Build buoypy for the NDBC stuff
WORKDIR /home/nru/repo/buoypy
RUN python setup.py install

# switch to the non-root user
USER nru

# change to the apsviz2 pipelines directory                        
WORKDIR /home/nru/repo/EDSASTAPPS/apsviz2

## set the python path
ENV PYTHONPATH=/home/nru/repo/AST:/home/nru/repo/EDSASTAPPS

##########
# The container is ready to accept the launch commands.
##########
