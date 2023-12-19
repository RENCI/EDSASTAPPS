#!/bin/bash -l
#
# run_harvester launcher script/wrapper. This script is intended to be used by our k8s system to generate 
# observation timeseries files for subsequent ingestion to the Harvester database. It is expected this script will
# be invoked daily by using cron.
#
# The git repo that downloads station and grid data muct be handled here
#
# Check environment

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
stoptime="$year-$month-$day $hour:00:00"

FINALDIR=$1

if [ -z "${LOG_PATH:-}" ]; then
   echo "LOG_PATH is not set: logs will be sent to $FINALDIR"
   LOG_PATH="$FINALDIR"
fi

if [ -z "${NDAYS:-}" ]; then
   echo "NDAYS is not set: Use default lookback behavior"
   NDAYSSET=""
else
    NDAYSSET="--ndays $NDAYS"
fi

# TEMP
#nru@4292a7f1007b:~/repo/EDSASTAPPS/run_harvester$ export HTTP_PROXY=http://proxy.renci.org:8080
#nru@4292a7f1007b:~/repo/EDSASTAPPS/run_harvester$ export HTTPS_PROXY=http://proxy.renci.org:8080
#nru@4292a7f1007b:~/repo/EDSASTAPPS/run_harvester$ export http_proxy=http://proxy.renci.org:8080
#nru@4292a7f1007b:~/repo/EDSASTAPPS/run_harvester$ export https_proxy=http://proxy.renci.org:8080

#
# git clone the grid data. The underlying directory structure is implied within the provided grid_to_stationfile_maps.yml file
#
git clone https://github.com/RENCI/AST_gridstations.git

#
# Where is Contrails authentication yml. Grab the secrets key from $CONTRAILS_KEY. Update the local secrets file
#
sed -i 's/xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx/'"$CONTRAILS_KEY"'/g' ./secrets/contrails.yml

#
# We only need to supply the proper URL to start the job
#

python ./run_fetch_pipeline_observation.py --stoptime "$stoptime"  --finalLOG "$LOG_PATH" --finalDIR "$FINALDIR"  --map_source_file './AST_gridstations/harvester_stations/sources_map.yaml' --contrails_auth './secrets/contrails.yml' $NDAYSSET

echo "Finished RUN_HARVESTER with status $?"
