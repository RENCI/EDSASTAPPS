#!/bin/bash
#SBATCH -t 128:00:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -J DailyHarvester
#SBATCH --mem-per-cpu 64000

## A quick script to process all preceding January days to now.


export PYTHONPATH=/projects/sequence_analysis/vol1/prediction_work/AST:/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS
export RUNTIMEDIR=/projects/ees/TDS/DataHarvesting/DAILY_HARVESTING
#export RUNTIMEDIR=./DEC-DAILIES

# Which PYTHON
PHOME="/projects/sequence_analysis/vol1/prediction_work/PythonMethods/anaconda3/bin"

# Where are the station files
STATIONDIR="/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/MAIN_GRID_DEFINITIONS"

# Where is Contrails authentication yml
CONT_AUTH="/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/run_harvester/secrets"

# Prepare invocation
cd /projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/run_harvester

# Observations

stime='2022-12-31 23:00:00'


# River data
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'CONTRAILS' --data_product 'river_water_level' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_rivers_largelist.csv"
# Contrails coastal data
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92 --data_source 'CONTRAILS' --data_product 'coastal_water_level' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_coastal.csv"
# Contrails River VOLUME
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'CONTRAILS' --data_product 'river_flow_volume' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_rivers.csv"
# Contrails air_pressure
# $PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'CONTRAILS' --data_product 'air_pressure' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_coastal.csv"

# NOAA coastal data
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NOAA' --data_product "water_level" --station_list "$STATIONDIR/noaa_stations.csv"
# NOAA predictions (tidal)
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NOAA' --data_product "predictions" --station_list "$STATIONDIR/noaa_stations.csv"
# NOAA Air pressure
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NOAA' --data_product "air_pressure" --station_list "$STATIONDIR/noaa_stations.csv"

# NDBC BUOY DATA
##$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NDBC' --data_product "wave_height" --station_list "$STATIONDIR/ndbc_buoys.csv"
# NDBC Wind speed
##$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NDBC' --data_product "wind_speed" --station_list "$STATIONDIR/ndbc_buoys.csv"
##$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NDBC' --data_product "air_pressure" --station_list "$STATIONDIR/ndbc_buoys.csv"


# NDBC_HISTORIC BUOY DATA
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NDBC_HISTORIC' --data_product "wave_height" --station_list "$STATIONDIR/ndbc_buoys.csv"
# NDBC Wind speed
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NDBC_HISTORIC' --data_product "wind_speed" --station_list "$STATIONDIR/ndbc_buoys.csv"
$PHOME/python run_get_observation.py --stoptime "$stime" --ndays -92  --data_source 'NDBC_HISTORIC' --data_product "air_pressure" --station_list "$STATIONDIR/ndbc_buoys.csv"
 
