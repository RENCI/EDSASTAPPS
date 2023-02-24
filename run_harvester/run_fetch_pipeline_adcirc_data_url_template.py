#!/usr/bin/env python
#
# A code is to facilitate the gathering of ADCIRC fort63/swan water data using the AST:harvester classes
# into a form for ingestion into our timeseries DB
#
# The expectation is this code is executed by a bash script within a cron 
#
# This wrapper only fetches ADCIRC data using the URL Template approach
# Only a fort63_style update is performed for this data execution
#
# No additional urls are attempted to be processed here

import os,sys
import shutil
import numpy as np
import pandas as pd
import datetime as dt
import harvester.get_adcirc_stations as get_adcirc_stations
import harvester.generate_urls_from_times as genurls # generate_urls_from_times
import harvester.fetch_adcirc_data as fetch_adcirc_data
import gridmap.grid_to_station_maps as grid_to_station_maps

from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities

varnamedict={}
varnamedict['fort.63']='zeta'
varnamedict['fort.61']='zeta'
varnamedict['swan_HS']='swan_HS'
varnamedict['swan_TPS']='swan_TPS'
varnamedict['swan_DIR']='swan_DIR'

def format_data_frames(df) -> pd.DataFrame:
    """
    A Common formatting used by all sources
    """
    df.index = df.index.strftime('%Y-%m-%dT%H:%M:%S')
    df.reset_index(inplace=True)
    df_out=pd.melt(df, id_vars=['TIME'])
    df_out.columns=('TIME','STATION',PRODUCT.upper())
    df_out.set_index('TIME',inplace=True)
    return df_out

def find_variable_name(url)->str:
    """
    """
    if "fort.61" in url: return varnamedict['fort.61']
    if "fort.63" in url: return varnamedict['fort.63']
    if "swan_HS" in url: return varnamedict['swan_HS']
    if "swan_TPS" in url: return varnamedict['swan_TPS']
    if "swan_DIR" in url: return varnamedict['swan_DIR']
    utilities.log.info(f'Failed to find varname {url}')
    sys.exit(1)

def change_ensemble_name(url, ensemble)->str:
    """
    """
    urlwords=url.split('/')
    urlwords[-2]=ensemble
    newurl='/'.join(urlwords)
    return newurl

def find_cast_status(df)->str:
    """
    Look into the provided metadata file to see if this was a declared NOWCAST or FORECAST
    (independant of the ensemble name)
    The relevant column name from fetch_stations is NAME=grid_type

    Parameter
        df: dataframe of station x NAME metadata
    Return
        cast_type: (str) either NOWCAST or FORECAST
    """
    data_list = df['NAME'].tolist()
    type_pair = list(set(data_list))[0] # fetchj_stations has already ensured these will all be the same grid and type
    utilities.log.info(type_pair)
    cast_type = type_pair.split('_')[-1]
    return cast_type.upper()

dformat='%Y-%m-%d %H:%M:%S'
GLOBAL_TIMEZONE='gmt' # Every source is set or presumed to return times in the zone
PRODUCT='water_level'

# Currently supported sources
SOURCES = ['ASGS']

##
## Run stations
##

def main(args):
    """
    Generally we anticipate inputting a STOPTIME
    Then the STARTTIME is ndays on the past
    """

    main_config = utilities.init_logging(subdir=None,config_file=os.path.join(os.path.dirname(__file__),'./config','main.yml'))

    # Grab the args
    map_file=args.map_file
    url = args.url
    return_sample_min=15
    ensemble=args.ensemble

    if url is None:
        utilities.log.error('For URL-based approaches template url must be specified')
        sys.exit(1)

    if map_file is None: 
        utilities.log.error('map_file must be specified')
        sys.exit(1)

    # Get the IO basic setep
    if args.finalDIR is None:
        rootdir=io_utilities.construct_base_rootdir(main_config['DEFAULT']['RDIR'], base_dir_extra=None)
    else:
        utilities.log.info('Override with finalDIR setting {}'.format(args.finalDIR))
        rootdir=io_utilities.construct_base_rootdir(args.finalDIR, base_dir_extra=None)

    if args.finalLOG is None:
        logdir=rootdir
    else:
        logdir=args.finalLOG

    utilities.log.info(f'Output logger directory specified to be {logdir}')

    if args.sources:
         utilities.log.info('Return list of sources')
         return SOURCES
         sys.exit(0)
    data_source = args.data_source
    data_product = args.data_product

    if data_source.upper() in SOURCES:
        utilities.log.info('Found selected data source {}'.format(data_source))
    else:
        utilities.log.error('Invalid data source {}'.format(data_source))
        sys.exit(1)

    url = args.url

    if args.ensemble is not None:
        url = change_ensemble_name(url, args.ensemble)
        utilities.log.info(f'Changing ensemble name to {args.ensemble}')

    variable_name = find_variable_name(url)
    utilities.log.info(f'Found variable name of {variable_name}')

    # Ensure the input URL is of the proper style (63 vs 61)
    if args.fort63_style:
        if variable_name is 'zeta':
            input_url=get_adcirc_stations.convert_urls_to_63style([url])[0] # Assumed to only have one of them
        else:
            input_url=get_adcirc_stations.convert_urls_to_63style_customfilename([url],filename='swan_HS.63.nc')[0] # Assumed to only have one of them
            utilities.log.info('SWAN style data file')
    else:
        if variable_name is 'zeta':
            input_url=get_adcirc_stations.convert_urls_to_61style([url])[0]
        else:
            utilities.log.info("No fort61_style approach for dara other than fort63/61.nc")
            sys.exit(1)
    utilities.log.info(input_url)
    urls=[input_url]

    # Fetch the grid specific station list.
    grid_name = fetch_adcirc_data.grab_gridname_from_url(urls)
    station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=grid_name, mapfile=map_file)
    utilities.log.info(f'Station file determined to be {station_file}')

    rpl = get_adcirc_stations.get_adcirc_stations(source='ASGS', product=args.data_product,
                station_list_file=station_file, fort63_style=args.fort63_style )

    data_adc,meta_adc=rpl.fetch_station_product(urls, return_sample_min=return_sample_min, fort63_style=args.fort63_style, variable_name=variable_name )

    # Get extra info for amending output filenames
    sitename=rpl.sitename
    gridname=rpl.gridname
    ensemble=rpl.ensemble
    starttime = dt.datetime.strptime(rpl.Tmin,'%Y%m%d%H').strftime("%Y-%m-%d %H:%M:%S")
    endtime = dt.datetime.strptime(rpl.Tmax,'%Y%m%d%H').strftime("%Y-%m-%d %H:%M:%S")
    utilities.log.info(f'Starttime {starttime}, endtime {endtime}')

    # Look inside meta and see what kind of data this is NOWCAST/FORECAST
    cast_type = find_cast_status(meta_adc)

    # Output the data files for subsequent ingesting
    df_adcirc_data = format_data_frames(data_adc)
    adcirc_metadata=sitename.upper()+'_'+ensemble.upper()+'_'+grid_name.upper()+'_'+cast_type+'_'+endtime.replace(' ','T')+'_'+starttime.replace(' ','T')+'_'+endtime.replace(' ','T')
    try:
        dataf=io_utilities.write_csv(df_adcirc_data, rootdir=rootdir,subdir='',fileroot='adcirc_stationdata',iometadata=adcirc_metadata)
        metaf=io_utilities.write_csv(meta_adc, rootdir=rootdir,subdir='',fileroot='adcirc_stationdata_meta',iometadata=adcirc_metadata)
        utilities.log.info('ADCIRC data has been stored {},{}'.format(dataf,metaf))
    except Exception as e:
        utilities.log.error('Error: ADCIRC: Failed Write {}'.format(e))
        sys.exit(1)

    log_time_meta = dt.datetime.strptime(rpl.Tmax,'%Y%m%d%H').strftime("%Y%m%d%H%M%S")
    utilities.log.info('Copy log file')
    shutil.copy(utilities.LogFile,'/'.join([logdir,f'run_harvester_adcirc_{log_time_meta}_log'])) 
    utilities.log.info('Finished with data source {}'.format(data_source))
    utilities.log.info('Finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--data_source', action='store', dest='data_source', default='ASGS', type=str,
                        help='choose supported data source: default = ASGS')
    parser.add_argument('--data_product', action='store', dest='data_product', default='water_level', type=str,
                        help='choose supported data product: default is water_level')
    parser.add_argument('--map_file', action='store', dest='map_file', default=None, type=str,
                        help='Location of the grid_to_stationfile_maps.yml data')
    parser.add_argument('--url', action='store', dest='url', default=None, type=str,
                        help='ASGS url to fetcb ADCIRC data')
    parser.add_argument('--ensemble', action='store',dest='ensemble', default=None,
                        help='str: Select appropriate ensemble Default is nowcast')
    parser.add_argument('--fort63_style', action='store_true', default=False,
                        help='Boolean: Will inform Harvester to use fort.63.methods to get station nodesids')
    parser.add_argument('--finalDIR', action='store', dest='finalDIR', default=None,
                        help='String: Custom location for the output dicts, PNGs and logs')
    parser.add_argument('--finalLOG', action='store', dest='finalLOG', default=None,
                        help='String: Custom location logs. If not specified logs go to the datadir')
    args = parser.parse_args()
    sys.exit(main(args))
