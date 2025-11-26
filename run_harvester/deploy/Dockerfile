##############
# stage 1: create a compressed conda virtual environment
##############
FROM continuumio/miniconda3 AS build

# author
LABEL author="Jeff Tilson"

# extra metadata
LABEL version="v0.0.0"
LABEL description="EDSASTAPPS Dockerfile."

# update conda
RUN conda update conda

# get the environment file
COPY ./environment.yml .

# Create the virtual environment
RUN conda env create -f environment.yml

# install conda pack to compress this stage
RUN conda install -c conda-forge conda-pack

# conpress the virtual environment
RUN conda-pack -n edsastapps -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# fix up the paths
RUN /venv/bin/conda-unpack

##############
# stage 2: create a python implementation using the stage 1 virtual environment
##############
FROM python:3.13-slim

# update the environment
RUN apt-get update

# install git
RUN apt-get -y install git

# clear out the apt cache
RUN apt-get clean

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# make the virtual environment active
ENV VIRTUAL_ENV=/venv
ENV PATH=/venv/bin:$PATH

# add user nru and switch to it
RUN useradd --create-home -u 1000 nru

# Copy in the rest of the code
RUN mkdir -p /home/nru/repo

# move the the code location
WORKDIR /home/nru/repo

# get the various source repos
RUN git clone https://github.com/RENCI/AST.git
RUN git clone https://github.com/RENCI/buoypy.git

# note this cloning will only get what has been pushed to GH:main and not anything in this local
# environment (including this file). suggest pushing this project to a branch, and specifying
# that branch in the git clone below.
RUN git clone https://github.com/RENCI/EDSASTAPPS.git

# move to the deploy directory
WORKDIR /home/nru/repo/EDSASTAPPS/run_harvester/deploy

# # make sure everything is read/write
RUN chmod 777 -R /home/nru

# switch to the buoypy directory and build it for the NDBC stations
WORKDIR /home/nru/repo/buoypy
RUN python setup.py install

# Now remove some issues that get dropped in
RUN pip uninstall -y protobuf

# switch to the non-root user
USER nru

# set the app working directory
WORKDIR /home/nru/repo/EDSASTAPPS/run_harvester

# testing only
# RUN git clone https://github.com/RENCI/AST_gridstations.git
# COPY main.yml .
# ENV URL="http://tds.renci.org/thredds/dodsC/2025/gfs/2025112200/NCSCv2.0/ht-ncfs.renci.org/NCSCv2.0_gfs/gfsforecast/fort.61.nc"
# ENV FINALDIR="./JUNK"
# ENV LOG_PATH="./JUNK"
# ENV RUNTIMEDIR=/tmp
# python run_fetch_pipeline_adcirc_data_multiprovider_url_template.py --data_source 'TDS' --url $URL --map_file './AST_gridstations/full_stationlist/grid_to_stationfile_maps.yml' --source_file './AST_gridstations/harvester_stations/sources_map.yaml'  --data_provider 'NOAA_STATIONS' --ensemble 'gfsforecast'
# python run_fetch_pegasus_observation_class.py --stoptime '2025-11-01 00:00:00' --map_source_file './sources_map.yaml' --contrails_auth './secrets/contrails.yml' --finalDIR ./JUNKJUNE --noaa_datum 'NAVD'

## set the python path
ENV PYTHONPATH=/home/nru/repo/AST:/home/nru/repo/EDSASTAPPS
