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
#SBATCH -J Reanalysis
#SBATCH --mem-per-cpu 64000

export PYTHONPATH=/projects/sequence_analysis/vol1/prediction_work/AST:/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS
export RUNTIMEDIR=./JUNK

YEAR=$1

GRID="hsofs"
DAILY=DAILY-$GRID-LP24

export INDIR=$RUNTIMEDIR/
export OUTROOT=$RUNTIMEDIR/$DAILY

##
## Test running the hsofs annual while including the knockout data: This defaults fort61_style
##

python compute_annual_errors.py --fort63_style --url "/projects/reanalysis/ADCIRC/ERA5/hsofs/$YEAR/fort.63.nc" --gridname 'hsofs' --map_file './supporting_data/grid_to_stationfile_maps.yml' --knockout  ./supporting_data/knockoutStation.json

##
## run the DAILY breakdown
##

python daily_lowpass_sampled_error.py  --inyear $YEAR  --input_directory $INDIR --output_directory $OUTROOT --yearly_file_properties 'yearly_file_properties.json' --main_yamlname './config/main.yml'

##
## run the interpolation 
##

python wrap_interpolation_daily_errorset.py --json_daily_file_errors $INDIR/daily_file_properties.json --input_directory $INDIR --output_directory $OUTROOT --main_yamlname './config/main.yml' --gridname $GRID --map_file './supporting_data/grid_to_stationfile_maps.yml'

