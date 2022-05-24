#
# Setup the basic env values relative to /projects/sequence_analysis/vol1/prediction_work/HARVESTOR/fetch_station_data
#
# This job is to be run every 6 hours to check for adcirc files
# Run the cron at 6 hours + 15 mins to help ensure the data actually got to the ASGS site
# We do not included $min:$sec values here

# NOTE: The [provided urls are as "templates" from which new urls are created

export PYTHONPATH=/projects/sequence_analysis/vol1/prediction_work/AST:/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS
#export RUNTIMEDIR=./JUNK
export RUNTIMEDIR=/projects/ees/TDS/DataHarvesting/DAILY_HARVESTING


# Which PYTHON
PHOME="/projects/sequence_analysis/vol1/prediction_work/PythonMethods/anaconda3/bin"

# Where is AST?
SRC="/projects/sequence_analysis/vol1/prediction_work/AST/run_harvester"

# Where are the station files
STATIONDIR="/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/MAIN_GRID_DEFINITIONS"

# Where is Contrails authentication yml
CONT_AUTH="/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/run_harvester/secrets"

# Prepare invocation
cd /projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/run_harvester

# Dig out the DAY, MONTH, and YEAR values. Always set hour:min:sec to 00:00:00

day="$(date +'%d')"
month="$(date +'%m')"
year="$(date +'%Y')"
hour="$(date +'%H')"
printf "Current rumtime in yyyy-mm-dd hh format %s\n" "$year-$month-$day $hour"

# Let's check to see if hour is either 00,06,12,or 18. Else do nothing and fix your cron launch parameters. *set to  

hour="00"
if [ "$hour" -eq "00" ] || [ "$hour" -eq "06" ] || [ "$hour" -eq "12" ] || [ "$hour" -eq "18" ]
    then
        # Get the nowcasts
        stoptime="$year-$month-$day $hour:00:00"
        $PHOME/python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.61.nc" --ensemble 'nowcast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays -2
        $PHOME/python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/ec95d/hatteras.renci.org/ec95d-nam-bob-da-nowcast/nowcast/fort.61.nc" --ensemble 'nowcast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays -2
        # Get the ec95d nowcasts
        $PHOME/python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/hsofs/hatteras.renci.org/hsofs-nam-bob-2021/nowcast/fort.61.nc" --ensemble 'namforecast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays -2
        #$PHOME/python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/ec95d/hatteras.renci.org/ec95d-nam-bob-da-nowcast/nowcast/fort.61.nc" --ensemble 'namforecast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays -2

        # NCSC_SAB_v1.23
        #$PHOME/python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/ncsc_sab_v1.23/hatteras.renci.org/ncsc123-nam-sb/nowcast/fort.61.nc" --ensemble 'nowcast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays -2
        #$PHOME/python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/ncsc_sab_v1.23/hatteras.renci.org/ncsc123-nam-sb/namforecast/fort.61.nc" --ensemble 'namforecast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays -2

        # Grab an ec05d from the PSC site as well
        #$PHOME/python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/ec95d/bridges2.psc.org/ec95d-nam-bob-psc/nowcast/fort.61.nc" --ensemble 'nowcast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays -2
        #$PHOME/python run_get_adcirc_data_url_template.py --stoptime "$stoptime" --data_source 'ASGS'  --url "http://tds.renci.org/thredds/dodsC/2022/nam/2022020200/ec95d/bridges2.psc.org/ec95d-nam-bob-psc/nowcast/fort.61.nc" --ensemble 'namforecast' --map_file './supporting_data/grid_to_stationfile_maps.yml' --ndays -2
       else
       # Cron should have taken care of this but here we are...
        echo "A non-6hourly time was specified. Cant do anything with this but exit"
fi
echo "Finished"

