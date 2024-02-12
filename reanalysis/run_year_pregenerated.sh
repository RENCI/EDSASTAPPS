#!/bin/sh

# 
# An example of how one might run a reanalysis pipeline. For this example we assume an input variable of YEAR. All output I/O for a given YEAR is amended such that
# if you run a caller script over many years, no clashes should occur.
# See the shell loopYears.sh for an example
#

#SBATCH -t 48:00:00
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1 
#SBATCH -p lowpri
#SBATCH -J Reanalysis
#SBATCH --mem-per-cpu 64000

#conda activate apsviz

export PYTHONPATH=/projects/prediction_work/AST:/projects/prediction_work/EDSASTAPPS

YEAR=$1

GRID="hsofs.V2"
DAILY=DAILY-$GRID-LP24

export RUNTIMEDIR=./$GRID/YEARLY-$YEAR
export INDIR=$RUNTIMEDIR/
export OUTROOT=$RUNTIMEDIR/$DAILY
export LOG_PATH=$RUNTIMEDIR

##
## Test running the hsofs annual while excluding knockout data: This defaults fort61_style
##

#inobsdir='/projects/prediction_work/NOAA.Reanalysis/43_OBSERVATIONS_BRIAN_DETRENDING'
#inobsfile='noaa_obs_1979_2022_detrended'

python compute_annual_errors_read_pregenerated_obs.py --fort63_style --url "/projects/reanalysis/ADCIRC/ERA5/$GRID/$YEAR/fort.63.d4.no-unlim.T.rc.nc" --gridname $GRID --map_file './supporting_data/grid_to_stationfile_maps.yml' --custom_fort63_name 'fort.63.d4.no-unlim.T.rc.nc' --obs_filename "/projects/prediction_work/NOAA.Reanalysis/43_OBSERVATIONS_BRIAN_DETRENDING/noaa_obs_1979_2022_detrended.csv" --meta_filename "/projects/prediction_work/NOAA.Reanalysis/43_OBSERVATIONS_BRIAN_DETRENDING/CERA_NOAA_HSOFS_stations_V4_withPR.csv"

##
## run the DAILY breakdown
##

##python daily_pregenerated_error.py  --inyear $YEAR  --input_directory $INDIR --output_directory $OUTROOT --yearly_file_properties 'yearly_file_properties.json' --main_yamlname './config/main.yml'

##
## run the interpolation 
##

##python wrap_interpolation_daily_errorset.py --json_daily_file_errors $INDIR/daily_file_properties.json --input_directory $INDIR --output_directory $OUTROOT --main_yamlname './config/main.yml' --gridname $GRID --map_file './supporting_data/grid_to_stationfile_maps.yml' --cv_testing

