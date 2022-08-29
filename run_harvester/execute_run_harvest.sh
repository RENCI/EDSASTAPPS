#!/bin/bash -eu
#
# run_harvester launcher script/wrapper. This script is intended to be used by our k8s system to generate 
# observation timeseries files for subsequent ingestion to the Harvester database. It is expected this script will
# be invoked daily by using cron.
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

day="$(date +'%d')"
month="$(date +'%m')"
year="$(date +'%Y')"
hour="$(date +'%H')"
printf "Current date in dd/mm/yyyy format %s\n" "$day-$month-$year $hour:00:00"
stoptime="$year-$month-$day 00:00:00"

FINALDIR=$1

#
# git clone the grid data. The underlying directory structure is implied within the provided grid_to_stationfile_maps.yml file

#git clone https://github.com/RENCI/AST_gridstations.git

#
# Where is Contrails authentication yml. Grab the secrets key from $CONTRAILS_KEY. Update the local secrets file
#

sed -i 's/xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx/'"$CONTRAILS_KEY"'/g' ./secrets/contrails.yml

#
# We only need to supply the proper URL to start the job
#

python ./run_fetch_pipeline_observation.py --stoptime "$stoptime"  --finalDIR "$FINALDIR"  --map_source_file './AST_gridstations/harvester_stations/sources_map.yaml' --contrails_auth './secrets/contrails.yml'

echo "Finished RUN_HARVESTER $URL with status $?"

