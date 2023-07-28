##############
# Docker file for the creation of the adcirc support tools (AST).
##############

FROM conda/miniconda3-centos7

# As it reads
RUN yum install -y git
RUN yum update
RUN yum install -y libxml2-devel
RUN yum install -y libxslt-devel
RUN yum install -y gcc
#RUN yum install google-chrome-stable.x86_64 --nogpgcheck

# get some credit
LABEL maintainer="jtilson@renci.org"

# # update conda and set the download channel. This will bost python to 3.10+ and update ssl certs code

RUN conda config --add channels conda-forge
RUN conda update -y conda

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
WORKDIR /home/nru/repo/EDSASTAPPS/run_harvester/deploy

# install the pip requirements
RUN pip install -r pip_requirements.txt

# # make sure everything is read/write
RUN chmod 777 -R /home/nru

# Build buoypy for the NDBC stuff
WORKDIR /home/nru/repo/buoypy
RUN python setup.py install

# switch to the non-root user
USER nru

WORKDIR /home/nru/repo/EDSASTAPPS/run_harvester

# This is just a reminder. SuperV will set these in k8s at launch time 
#export http_proxy=http://proxy.renci.org:8080
#export https_proxy=http://proxy.renci.org:8080
#export HTTP_PROXY=http://proxy.renci.org:8080
#export HTTPS_PROXY=http://proxy.renci.org:8080

## set the python path
ENV PYTHONPATH=/home/nru/repo/AST:/home/nru/repo/EDSASTAPPS

##########
# The container is ready to accept the launch commands.
##########
