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

if [ -z "${RUNTIMEDIR:-}" ]; then
   echo "RUNTIMEDIR is not set"
fi

if [ -z "${CONTRAILS_KEY:-}" ]; then
   echo "CONTRAILS_KEY is not set: Abort"
   exit 1
fi

# Check input varibles

URL=$1
GRID=$2
FINALDIR=$3
INSTANCE=$4

#
# git clone the grid data. The underlying directory structure is implied within the provided grid_to_stationfile_maps.yml file

git clone https://github.com/RENCI/AST_gridstations.git

#
# Where is Contrails authentication yml. Grab the secrets key from $CONTRAILS_KEY. Update the local secrets file
#

sed -i 's/xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx/'"$CONTRAILS_KEY"'/g' ./secrets/contrails.yml

#
# We only need to supply the proper URL to start the job
#

python ./apsviz_station_products.py --fort63_style --construct_jsons --url "$URL" --ndays -4  --return_sample_min 60 --gridname "$GRID" --ensemble='nowcast'  --map_file './AST_gridstations/full_stationlist/grid_to_stationfile_maps.yml' --contrails_auth ./secrets/contrails.yml --finalDIR "$FINALDIR" --instanceId "$INSTANCE" |& tee "$INSTANCE"_stdout

echo "Finished OBSMOD $URL with status $?"

