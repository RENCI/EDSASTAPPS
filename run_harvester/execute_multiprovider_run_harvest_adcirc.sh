#/bin/bash -l
#
# launcher adcirc harvester script. This script is intended to be used by our k8s apsviz2 pipeline
#
# The input arguments have been expanded in this new launcher. Setting finalDir, will override any value of $RUNTIMEDIR
# The log file will be copied to finalDir at the end of the processing
#
# The git repo that downloads station and grid data must be handled here
#
# The optional INSTANCEID is a presumptive unique runtime tag that gets added to the output data and metadata filenames
# if INSTANCEID is not set (Null) then it is quietly ignored
#
# Check environment

URL=$1

if [ -z "${PYTHONPATH:-}" ]; then
   echo "PYTHONPATH is not set"
   exit 1
fi

FINALDIR=$2

if [ -z "${LOG_PATH:-}" ]; then
   echo "LOG_PATH is not set: logs will be sent to $FINALDIR"
   LOG_PATH="$FINALDIR"
fi

PROVIDER=$3

INSTANCEID=$4
if [ -z "$INSTANCEID" ]
then
    echo "No INSTANCEID"
else
    echo "Got an INSTANCEID of $INSTANCEID"
fi


# git clone the grid data. The underlying directory structure is implied within the provided grid_to_stationfile_maps.yml file
git clone https://github.com/RENCI/AST_gridstations.git


# Run the actual input URL first without changing the ensemble name

echo "Begin ADCIRC Forecast fetch"
python run_fetch_pipeline_adcirc_data_multiprovider_url_template.py --data_source 'TDS'  --url "$URL" --map_file  './AST_gridstations/full_stationlist/grid_to_stationfile_maps.yml'  --source_file './AST_gridstations/harvester_stations/sources_map.yaml' --data_provider "$PROVIDER" --finalDIR "$FINALDIR" --finalLOG "$LOG_PATH" --modelid "$INSTANCEID"
echo "Finished Forecast ADCIRC fetch"

# Now we want to rerun the associated nowcast for thie input URL

echo "Begin ADCIRC Nowcast fetch"
python run_fetch_pipeline_adcirc_data_multiprovider_url_template.py --data_source 'TDS'  --url "$URL" --map_file  './AST_gridstations/full_stationlist/grid_to_stationfile_maps.yml'  --source_file './AST_gridstations/harvester_stations/sources_map.yaml' --data_provider "$PROVIDER" --finalDIR "$FINALDIR" --finalLOG "$LOG_PATH" --ensemble "nowcast" --modelid "$INSTANCEID"
echo "Finish ADCIRC Nowcast fetch"

echo "Finished FETCH ADCIRC $URL with status $?"

