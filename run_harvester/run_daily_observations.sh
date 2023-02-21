#!/bin/bash
#SBATCH -t 512:00:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -J DailyHarvester
#SBATCH --mem-per-cpu 64000

## A quick script to process all preceding January days to now.

export PYTHONPATH=/projects/sequence_analysis/vol1/prediction_work/AST:/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS
export RUNTIMEDIR=/projects/ees/TDS/DataHarvesting/DAILY_HARVESTING
#export RUNTIMEDIR=./DAILIES

# Which PYTHON
PHOME="/projects/sequence_analysis/vol1/prediction_work/PythonMethods/anaconda3/bin"

# Where are the station files
STATIONDIR="/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/MAIN_GRID_DEFINITIONS"

# Where is Contrails authentication yml
CONT_AUTH="/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/run_harvester/secrets"

# Prepare invocation
cd /projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/run_harvester

# Observations

#for DAYS in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26; do
for MONTHS in 01 02 03 04 05 06 07 08 09 10 11 12; do
    #for DAYS in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 31; do
    for DAYS in 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 31; do
        # River data
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'CONTRAILS' --data_product 'river_water_level' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_rivers.csv"
        # Contrails coastal data
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00' --data_source 'CONTRAILS' --data_product 'coastal_water_level' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_coastal.csv"
        # Contrails River VOLUME
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'CONTRAILS' --data_product 'river_flow_volume' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_rivers.csv"
        # Contrails air_pressure
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'CONTRAILS' --data_product 'air_pressure' --config_name "$CONT_AUTH/contrails.yml" --station_list "$STATIONDIR/contrails_stations_coastal.csv"

        # NOAA coastal data
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NOAA' --data_product "water_level" --station_list "$STATIONDIR/noaa_stations.csv"
        # NOAA predictions (tidal)
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NOAA' --data_product "predictions" --station_list "$STATIONDIR/noaa_stations.csv"
        # NOAA Air pressure
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NOAA' --data_product "air_pressure" --station_list "$STATIONDIR/noaa_stations.csv"

        # NDBC BUOY DATA
        #$PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NDBC' --data_product "wave_height" --station_list "$STATIONDIR/ndbc_buoys.csv"
        # NDBC Wind speed
        #$PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NDBC' --data_product "wind_speed" --station_list "$STATIONDIR/ndbc_buoys.csv"
        #$PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NDBC' --data_product "air_pressure" --station_list "$STATIONDIR/ndbc_buoys.csv"


        # NDBC_HISTORIC BUOY DATA
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NDBC_HISTORIC' --data_product "wave_height" --station_list "$STATIONDIR/ndbc_buoys.csv"
        # NDBC Wind speed
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NDBC_HISTORIC' --data_product "wind_speed" --station_list "$STATIONDIR/ndbc_buoys.csv"
        $PHOME/python run_get_observation.py --stoptime '2022-'$MONTHS'-'$DAYS' 00:00:00'  --data_source 'NDBC_HISTORIC' --data_product "air_pressure" --station_list "$STATIONDIR/ndbc_buoys.csv"




 
    done
done
