#!/bin/bash
#SBATCH -t 128:00:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -J ADCIRCHarvester
#SBATCH --mem-per-cpu 64000

## A quick script to process all preceding January days to now.

export PYTHONPATH=/projects/sequence_analysis/vol1/prediction_work/AST:/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS
export RUNTIMEDIR=/projects/ees/TDS/DataHarvesting/DAILY_HARVESTING

# ADCIRC

for MONTHS in 04; do
    for DAYS in  01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 31; do
        for HOURS in 00 06 12 18 ; do
            stoptime="2022-$MONTHS-$DAYS $HOURS:00:00" 
            echo "$stoptime"
            python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.61.nc" --ensemble 'nowcast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays 0 
            python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.61.nc" --ensemble 'namforecast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays 0 
            #python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/ec95d/hatteras.renci.org/ec95d-nam-bob-da-nowcast/nowcast/fort.61.nc" --ensemble 'nowcast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays 0 
         done
    done
done

