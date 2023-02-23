#!/usr/bin/env python
# This is a followon to the original run_get_observation.py code. This version is much more amenable to handling multiple data_products.
# This version, reads a yasml from the AST_gridnames repo to determine which sources to fetch and which data_products to fetch from each 
# source
#
# We intentionally have time ranges that are overlapping
#
# The filenames are going to have the following nomenclature:
#
# ensemble: An arbitrary string. According to a cursory look at tds, this has values such as:
#    nowcast, nhc0fcl, veerright, etc. So we will set the following defaults:
# grid: hsofs,ec95d,etc
#
#


#
# Here we are simulating having run a series of fetches from the Harvester and storing the data to csv files.
# These files will be used for building and testing out schema for the newq DB 
#
# We intentionally have time ranges that are overlapping
#
# The filenames are going to have the following nomenclature:
#
# ensemble: An arbitrary string. According to a cursory look at tds, this has values such as:
#    nowcast, nhc0fcl, veerright, etc. So we will set the following defaults:
# grid: hsofs,ec95d,etc
#
#

import os,sys
import shutil
import numpy as np
import pandas as pd
import datetime as dt
import harvester.get_observations_stations as get_obs_stations
from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities

#$PHOME/python  "$SRC"/fetch_data.py --data_source 'NOAA' --data_product 'water_level' --stoptime "$stoptime" --station_list "$STATIONDIR"/noaa_stations.csv --ofile "$RUNTIMEDIR" --ometafile "$RUNTIMEDIR" 

def format_data_frames(df,product) -> pd.DataFrame:
    """
    A Common formatting used by all sources
    """
    prod_name = PRODUCT[product].upper() if product in PRODUCT.keys() else 'NOAME'
    print(f'Data Product: {PRODUCT[product]}')
    df.index = df.index.strftime('%Y-%m-%dT%H:%M:%S')
    df.reset_index(inplace=True)
    df_out=pd.melt(df, id_vars=['TIME'])
    df_out.columns=('TIME','STATION',prod_name)
    df_out.set_index('TIME',inplace=True)
    return df_out

##
## Globals
##

dformat='%Y-%m-%d %H:%M:%S'
GLOBAL_TIMEZONE='gmt' # Every source is set or presumed to return times in the zone
##PRODUCT='water_level'

PRODUCT={
        'water_level':'water_level',
        'hourly_height':'water_level',
        'predictions':'water_level',
        'river_water_level':'water_level',
        'river_flow_volume':'river_flow_volume',
        'coastal_water_level':'water_level',
        'wave_height':'water_level',
        'air_pressure':'air_pressure',
        'wind_speed':'wind_speed'
         }

# Currently supported sources
SOURCES = ['NOAA','CONTRAILS','NDBC']

##
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

##
## Select where to write these files
##
    if args.finalDIR is None:
        rootdir=io_utilities.construct_base_rootdir(config['DEFAULT']['RDIR'], base_dir_extra=None)
    else:
        print('Override with finalDIR setting {}'.format(args.finalDIR))
        rootdir=io_utilities.construct_base_rootdir(args.finalDIR, base_dir_extra=None)

    if args.finalLOG is None:
        logdir=rootdir
    else:
        logdir=args.finalLOG

    utilities.log.info(f'Output directory specified to be {rootdir}')
    utilities.log.info(f'Output logger directory specified to be {logdir}')

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

##
## Build stoptime metadata for output logfle annotation
##
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

    # Process all products for the given source
    for data_source in data_sources:
        output_fileroot,output_metafileroot = construct_root_filenames(data_source)
        source_products = source_config[data_source]['SOURCES'].values() 
        station_file=source_config[data_source]['STATION_FILE']
        data_source_short = 'CONTRAILS' if 'CONTRAILS' in data_source else data_source
        output_fileroot,output_metafileroot = construct_root_filenames(data_source_short)
        utilities.log.info(f'File rootnames are {output_fileroot} and {output_metafileroot}')
        for data_product in source_products:
            print(f'Processing source {data_source} and product {data_product}')
            metadata = construct_metadata(data_product, endt)
            utilities.log.info(f'Preparing to fetch {data_product} from {data_source}')
            obs = get_obs_stations.get_obs_stations(source=data_source_short.upper(), product=data_product,
                contrails_yamlname=args.contrails_auth, knockout_dict=None, station_list_file=station_file)
            # Get data at highest resolution. Return at 15min intervals
            try: # DO this in case we got no data from any station
                data,meta=obs.fetch_station_product(time_range, return_sample_min=15, interval='None' )
                df_data = format_data_frames(data, data_product) # Melt data to support later database updates
                dataf,metaf = write_data_todisk(rootdir,data_source_short,df_data,meta,metadata,output_fileroot,output_metafileroot)
            except IndexError as e:
                utilities.log.info(f'No data products {data_product} for source {data_source}. Skip')
                pass
            utilities.log.info(f'Finished with data source {data_source}')

##
## Copy over logfile into the rootdir so as to reside with the output data sets
##

    shutil.copy(utilities.LogFile,'/'.join([logdir,f'run_harvester_obs_{log_time_meta}_log'])) # Copy and rename to logs for apsviz2 pipeline to find
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
    parser.add_argument('--map_source_file', action='store',dest='map_source_file', default=None,
                        help='str: Select appropriate map_source_file yml for source processing list')
    parser.add_argument('--finalDIR', action='store', dest='finalDIR', default=None,
                        help='String: Custom location for the output dicts, PNGs and (potentially) logs')
    parser.add_argument('--finalLOG', action='store', dest='finalLOG', default=None,
                        help='String: Custom location logs. If not specified logs go to the datadir')
    args = parser.parse_args()
    sys.exit(main(args))
