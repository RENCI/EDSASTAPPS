#
# Setup the basic env values relative to /projects/prediction_work/HARVESTOR/fetch_station_data

export PYTHONPATH=/projects/prediction_work/AST:/projects/prediction_work/EDSASTAPPS
export RUNTIMEDIR=./DAILY-TEST-NEW
#export RUNTIMEDIR=/projects/ees/TDS/DataHarvesting/DAILY_HARVESTING

# Which PYTHON
PHOME="/projects/prediction_work/PythonMethods/anaconda3/bin"

# Where is AST?
SRC="/projects/prediction_work/AST/harvester"

# Where are the station files
STATIONDIR="/projects/prediction_work/EDSASTAPPS/MAIN_GRID_DEFINITIONS"

# Where is Contrails authentication yml
CONT_AUTH="/projects/prediction_work/EDSASTAPPS/run_harvester/secrets"

# Where is my WD 
cd /projects/prediction_work/EDSASTAPPS/run_harvester

day="$(date +'%d')"
month="$(date +'%m')"
year="$(date +'%Y')"
hour="$(date +'%H')"
printf "Current date in dd/mm/yyyy format %s\n" "$day-$month-$year $hour:00:00"
stoptime="$year-$month-$day 00:00:00"

# River data
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'CONTRAILS' --data_product 'river_water_level' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_rivers.csv"

# Contrails coastal data
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'CONTRAILS' --data_product 'coastal_water_level' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_coastal.csv"

# Contrails River Volume
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'CONTRAILS' --data_product 'river_flow_volume' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_rivers.csv"

# Contrails River Stream Elevation
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'CONTRAILS' --data_product 'river_stream_elevation' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_rivers.csv"

# NOAA coastal data
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'NOAAWEB' --data_product "water_level" --station_list "$STATIONDIR/noaa_stations.csv"

# NOAA predictions (tidal)
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'NOAAWEB' --data_product "predictions" --station_list "$STATIONDIR/noaa_stations.csv"

# NDBC buoy data
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'NDBC' --data_product "wave_height" --station_list "$STATIONDIR/ndbc_buoys.csv"

# NOAA air pressure
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'NOAAWEB' --data_product "air_pressure" --station_list "$STATIONDIR/noaa_stations.csv"

# NDBC wind speed 
$PHOME/python run_get_observation.py --stoptime "$stoptime" --data_source 'NDBC' --data_product "wind_speed" --station_list "$STATIONDIR/ndbc_buoys.csv"


echo "Finished"
