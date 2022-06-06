##############
# Docker file for the creation of the adcirc support tools (AST).
##############
FROM continuumio/miniconda3

# get some credit
LABEL maintainer="jtilson@renci.org"

# make sure the container is up-to-date
RUN apt-get update

# Start off getting Python 3.7.3
RUN conda install -y python=3.7.3

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
#WORKDIR /repo

# get the AST repo - needs to be made public first
RUN git clone https://github.com/RENCI/AST.git
RUN git clone https://github.com/RENCI/EDSASTAPPS.git
RUN git clone https://github.com/nickc1/buoypy.git

# move to the deploy directory
WORKDIR /home/nru/repo/EDSASTAPPS/apsviz2/deploy
#WORKDIR /repo/EDSASTAPPS/apsviz2/deploy
                                
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
#WORKDIR /repo/EDSASTAPPS/apsviz2

## set the python path
ENV PYTHONPATH=/home/nru/repo/AST:/home/nru/repo/EDSASTAPPS:/home/nru/repo/buoypy
#ENV PYTHONPATH=/repo/AST:/repo/EDSASTAPPS

# set the location of the output directory
#ENV RUNTIMEDIR=/data

##########
# The container is ready to accept the launch commands.
##########
