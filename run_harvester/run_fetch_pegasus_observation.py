#!/usr/bin/env python

# This is a modification of the original run_get_observation.py code better suited for execution within a Pegasus workflow. 
# This version is passed a single source yaml file whch points to a singe station list file. Both files must be local to the invocation directory.
# Though the yaml continues to support specifying multipe data products, to increase parallelism, it is suggested to have one data product per yaml.
# Sample yaml/station lists files may be downloaded from the RENCI AST_gridnames github repo and looking in the observations directory
#
# A user may specify at the CLI to fetch either the raw data only or return daa that has been smoothed using a rolling average. The output filename
# metadata change slomewhat to reflect this choice.
#
# If the user chooses to return SMOOTHED data the following procedure occurs:
#    a) The Raw data is fetched at best sampling_min depending on the data
#    b) A centered rolling average is applied.
#    c) The averaged data interpolated as a linear function to fill only single remining missing gaps. The endpoints.
#    d) The data are resampled to the input sampling_min

#
# OUTPUT products
# There are two primary OUTPUT files generated
# Product data. These are (Time x station) dataframes containing the value of the product in question
# Meta data.    These report characteristcs of each station. Such as names, longitude, latitude, etc.
# Completely empty stations over the chosen time range will be removed from the data set

# The output filenames will have the following nomenclature
# Product data: <source>_stationdata_<product name>_<time mark>{_smoothed}.csv
# Meta data: <source>_stationdata_meta_<product name>_<time mark>{_smoothed}.csv
#
# Examples for a NOAAWEB invocation to fetch water_level, starting from 2025-08-01T00:00:00 and looking back two days (--ndays --2) will be
#    noaaweb_stationdata_water_level_2025-08-01T00:00:00.csv
#    noaaweb_stationdata_meta_water_level_2025-08-01T00:00:00.csv
#
# Smoothed version would be 
#    noaaweb_stationdata_water_level_2025-08-01T00:00:00_smoothed.csv
#    noaaweb_stationdata_meta_water_level_2025-08-01T00:00:00_smoothed.csv

# A typical invocation to fetch (default) 2 days worth of NOAA station data beginning at 2025-08-01 00:00:00, smoothing with a rolling average window of 13 units
# while fetching raw data on a 60 min sampling. So the rolling average will conssst of a window 11 (60min steps) wide. Output files will be placed into finalDIR.
# might look as follows

# python run_fetch_pegasus_observation.py --stoptime "2025-08-01 00:00:00"  --map_source_file './NOAA_sources.yaml' 
#      --contrails_auth './secrets/contrails.yml' --finalDIR "./output" --smooth --nwindow 11 --sampling_min 60

import os,sys
import shutil
import numpy as np
import pandas as pd
import datetime as dt
import harvester.get_observations_stations as get_obs_stations
from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities

## Globals

dformat='%Y-%m-%d %H:%M:%S'

PRODUCT={
        'water_level':'water_level',
        'hourly_height':'water_level',
        'predictions':'water_level',
        'river_water_level':'water_level',
        'river_stream_elevation':'water_level',
        'river_flow_volume':'river_flow_volume',
        'coastal_water_level':'water_level',
        'wave_height':'water_level',
        'air_pressure':'air_pressure',
        'wind_speed':'wind_speed'
         }

SOURCES = ['USGS','USGS_RIVERS','NOAA','CONTRAILS','NDBC','NDBC_HISTORIC']
DICT_SOURCES={'USGS_RIVERS':'USGS_RIVERS','USGS':'USGS','NOAA_STATIONS':'NOAAWEB','NDBC_HISTORIC':'NDBC_HISTORIC','NDBC_BUOYS':'NDBC','CONTRAILS_COASTAL':'CONTRAILS','CONTRAILS_RIVERS':'CONTRAILS'}

def construct_root_filenames(data_source):
    """
    """
    output_fileroot=f'{data_source.lower()}_stationdata'
    output_metafileroot=f'{data_source.lower()}_stationdata_meta'
    return output_fileroot,output_metafileroot

def construct_metadata(data_product, endt):
    """   
    A convenience function used by the Harvester ingest processor
    """
    if 'river_water_level' in data_product:
         metadata=f'{data_product}_{endt}_RIVERS'
    elif 'river_stream_elevation' in data_product:
         metadata=f'{data_product}_{endt}_RIVERS'
    elif 'coastal_water_level' in data_product:
         metadata=f'{data_product}_{endt}_COASTAL'
    else:
        metadata=f'{data_product}_{endt}'
    return metadata

def write_data_todisk(rootdir,data_source,data,meta,metadata,output_fileroot,output_metafileroot): 
    """
    """
    # Output
    try:
        dataf=io_utilities.write_csv(data, rootdir=rootdir,subdir='',fileroot=output_fileroot,iometadata=metadata)
        utilities.log.info(f'{data_source.upper()} data has been stored {dataf}')
        metaf=io_utilities.write_csv(meta, rootdir=rootdir,subdir='',fileroot=output_metafileroot,iometadata=metadata)
        utilities.log.info(f'{data_source.upper()} metadata has been stored {metaf}')
    except Exception as e:
        utilities.log.error('Error:{} Failed Write {}'.format(data_source,e))
        sys.exit(1)
    return dataf,metaf

##
## Run stations
##

def main(args):
    """
    Generally we anticipate inputting a STOPTIME
    Then the STARTTIME is ndays on the past
    """

    # Set logging
    main_config = utilities.init_logging(subdir=None,config_file=os.path.join(os.path.dirname(__file__),'./config','main.yml'))
    utilities.log.info("Product Level Working in {}.".format(os.getcwd()))

## Select where to write these files

    rootdir=io_utilities.construct_base_rootdir(args.finalDIR, base_dir_extra=None)
    if args.finalLOG is None:
        logdir=rootdir
    else:
        logdir=args.finalLOG

    sampling_min=args.sampling_min

    ## If choosing smooth then only smoothed data are returned
    do_smooth = args.smooth
    if do_smooth:
        nwindow=args.nwindow
        utilities.log.info(f'Smoothing window requested {nwindow}')

    utilities.log.info(f'Requested sampling frequency {sampling_min}')
    utilities.log.info(f'Output directory specified to be {rootdir}')
    utilities.log.info(f'Output logger directory specified to be {logdir}')
    utilities.log.info(f'Output user directive for smoothing {do_smooth}')

    # Get list of sources/products to process
    source_config = utilities.load_config(args.map_source_file)['SOURCEMAP']

    # Setup times range tuple
    if args.stoptime is not None:
        time_stop=dt.datetime.strptime(args.stoptime,dformat)
    else:
        time_stop=dt.datetime.now()
    time_start=time_stop+dt.timedelta(days=args.ndays) # How many days BACK
    starttime=dt.datetime.strftime(time_start, dformat)
    endtime=dt.datetime.strftime(time_stop, dformat)
    utilities.log.info('Selected time range is {} to {}, ndays is {}'.format(starttime,endtime,args.ndays))

## Build stoptime metadata for output logfle annotation
    log_time_meta = dt.datetime.strftime(time_stop,'%Y%m%d%H%M%S')

##
## Copy log file to finalDir
##
    shutil.copy(utilities.LogFile,'/'.join([rootdir,f'logs'])) # Copy and rename to logs for apsviz2 pipeline to find
    utilities.log.info('Copy log file')

    time_range=(starttime,endtime)

    # metadata are used to augment filename
    endt=endtime.replace(' ','T')

    # Process the sources/products from the map_source_file
    data_sources = list(source_config.keys())

    raw_sampling_min = sampling_min if not do_smooth else 0

    # Process all products for the given source
    for data_source_in in data_sources:
        data_source=DICT_SOURCES[data_source_in]
        utilities.log.info(f'Found the provider name for {data_source}')

        output_fileroot,output_metafileroot = construct_root_filenames(data_source_in)
        source_products = source_config[data_source_in]['SOURCES'].values() 

        station_file=source_config[data_source_in]['STATION_FILE']
        print(f' station file {station_file}')

        data_source_short = 'CONTRAILS' if 'CONTRAILS' in data_source_in else DICT_SOURCES[data_source_in]
        output_fileroot,output_metafileroot = construct_root_filenames(data_source_short)
        utilities.log.info(f'File rootnames are {output_fileroot} and {output_metafileroot}')
        for data_product in source_products:
            if data_product=='predictions':
                print(f' {data_product} type. Need to extend stop time by 120 hours')
                ext_time_stop=time_stop+dt.timedelta(hours=120)
                ext_starttime=dt.datetime.strftime(time_start, dformat)
                ext_endtime=dt.datetime.strftime(ext_time_stop, dformat)
                time_range_use=(ext_starttime,ext_endtime)
                utilities.log.info('Force extended time range to {} and {}'.format(ext_starttime,ext_endtime))
            else:
                time_range_use=time_range
            metadata = construct_metadata(data_product, endt)
            utilities.log.info(f'Preparing to fetch {data_product} from {data_source}')

            ## Requires RAW data even if only wanting the smoothed data
            obs = get_obs_stations.get_obs_stations(source=data_source_short.upper(), product=data_product,
                contrails_yamlname=args.contrails_auth, knockout_dict=None, station_list_file=station_file)

            ## Expose option to SMOOTH data to the desired window width
            ## If smooth is selected choose the full resolution to get raw data, smooth and then reset freq
            try: # Do this in case we got no data from any station
                data,meta=obs.fetch_station_product(time_range_use, return_sample_min=raw_sampling_min)
                df_data = data # format_data_frames(data, data_product) # Melt data to support later database updates
                ##df_data=df_data.set_index('TIME')
                df_data.index = pd.to_datetime(df_data.index) # Do this now and it propogates to the subsequent daasets
                print(f' df RAW {raw_sampling_min}: {df_data}')
                if do_smooth:
                    data_thresholded = obs.remove_missingness_stations(data, max_nan_percentage_cutoff=100) # (maxmum percentage of allowable nans)
                    meta_list = set(data_thresholded.columns.tolist()).intersection(meta.index.to_list())
                    meta_thresholded = meta.loc[meta_list]
                    data_smoothed = obs.fetch_smoothed_station_product(data_thresholded, return_sample_min=sampling_min, window=nwindow)
                    print('Got the smoothed')
                    print(data_smoothed)
                    print('finished the print smoothed')
                    df_data_smooth = data_smoothed # format_data_frames(data_smoothed, data_product) 
                    dataf,metaf = write_data_todisk(rootdir,data_source_short,df_data_smooth,meta_thresholded,f'{metadata}_smoothed',output_fileroot,output_metafileroot)
                else:
                    dataf,metaf = write_data_todisk(rootdir,data_source_short,df_data,meta,metadata,output_fileroot,output_metafileroot)
            except Exception as e: # IndexError as e:
                utilities.log.info(f'No data products {data_product} for source {data_source}. Skip')
                pass
            utilities.log.info(f'Finished with data source {data_source}')

##
## Copy over logfile into the rootdir so as to reside with the output data sets
##

    shutil.copy(utilities.LogFile,'/'.join([logdir,f'run_harvester_obs_{log_time_meta}.log'])) # Copy and rename to logs for apsviz2 pipeline to find
    utilities.log.info('Copy log file')

    utilities.log.info('Finished') 

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--ndays', action='store', dest='ndays', default=-2, type=int,
                        help='Number of look-back days from stoptime (or now): default -1')
    parser.add_argument('--stoptime', action='store', dest='stoptime', default=None, type=str,
                        help='Desired stoptime YYYY-mm-dd HH:MM:SS. Default=now')
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--contrails_auth', action='store', dest='contrails_auth', default=None, type=str,
                        help='Choose a non-default contrails auth contrails_auth')
    parser.add_argument('--map_source_file', required=True, action='store',dest='map_source_file', default=None,
                        help='str: Select appropriate map_source_file yml for source processing list')
    parser.add_argument('--finalDIR', required=True, action='store', dest='finalDIR', default=None,
                        help='String: Custom location for the output dicts, PNGs and (potentially) logs')
    parser.add_argument('--finalLOG', action='store', dest='finalLOG', default=None,
                        help='String: Custom location logs. If not specified logs go to the datadir')
    parser.add_argument('--sampling_min', action='store', dest='sampling_min', default=15, type=int,
                        help='Returned data sampling frequency: default 15 (min)')
    parser.add_argument('--smooth', action='store_true', default=False,
                        help='Boolean: Will inform Harvester to return smoothed data')
    parser.add_argument('--nwindow', action='store', dest='nwindow', default=11, type=int,
                        help='Rolling window width in steps of sapling_min units. If smooth is selecte: default 11')

    args = parser.parse_args()
    sys.exit(main(args))
