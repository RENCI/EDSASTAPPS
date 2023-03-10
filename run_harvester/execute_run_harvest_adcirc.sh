#!/bin/bash -eu
#
# apsviz launcher script/wrapper. This script is intended to be used by our k8s apsviz2 pipeline to generate 
# station/buoy insets figures.
#
# The inpout argumenbts have become much simpler in this new launcher. setting finalDir, will override any value of $RUNTIMEDIR
# The log file will be copied to finalDir at the end of the processing
#
# The git repo that downloads station and grid data muct be handled here
#
# Check environment

if [ -z "${PYTHONPATH:-}" ]; then
   echo "PYTHONPATH is not set"
   exit 1
fi

FINALDIR=$2

if [ -z "${LOG_PATH:-}" ]; then
   echo "LOG_PATH is not set: logs will be sent to $FINALDIR"
   LOG_PATH="$FINALDIR"
fi

#if [ -z "${NDAYS:-}" ]; then
#   echo "NDAYS is not set: Use default lookback behavior"
#   NDAYSSET=""
#else
#   echo "Override input NDAYSSET. Reset to zero"
#   NDAYSSET="--ndays 0"
#fi

URL=$1

#
# git clone the grid data. The underlying directory structure is implied within the provided grid_to_stationfile_maps.yml file

git clone https://github.com/RENCI/AST_gridstations.git


# Need to run the actual input URL first without changing the ensemble name

python run_fetch_pipeline_adcirc_data_url_template.py --data_source 'ASGS'  --url "$URL" --map_file './AST_gridstations/full_stationlist/grid_to_stationfile_maps.yml' --fort63_style --finalDIR "$FINALDIR" --finalLOG "$LOG_PATH" 

# Now we want to rerun the associated nowcast for thie input URL

python run_fetch_pipeline_adcirc_data_url_template.py --data_source 'ASGS'  --url "$URL" --ensemble='nowcast' --map_file './AST_gridstations/full_stationlist/grid_to_stationfile_maps.yml' --fort63_style  --finalLOG "$LOG_PATH"  --finalDIR "$FINALDIR" 

echo "Finished ADCIRC $URL with status $?"

