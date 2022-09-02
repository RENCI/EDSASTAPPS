#!/usr/bin/env python
#
# A small code to facilitate the gathering of ADCIRC water data using the AST:harvester classes
# The expectation is this code is executed by a bash script within a cron 
#
# This wrapper only fetches ADCIRC data using the URL Template approach
# Only a fort61_style update is performed for this data execution
#
# This approach expects to fetch urlsa based on the url_framework.yml

#python get_adcirc_stations.py --url "http://tds.renci.org/thredds/dodsC/2021/al09/11/hsofs/hatteras.renci.org/hsofs-al09-bob/nhcOfcl/fort.61.nc" --data_source 'ASGS' --timeout 8  --fort63_style


## NEED to detemrine what the variable_name is ad pass it to fetch_adcirc_stations()

import os,sys
import numpy as np
import pandas as pd
import datetime as dt
import harvester.get_adcirc_stations as get_adcirc_stations
import harvester.generate_urls_from_times as genurls # generate_urls_from_times
import harvester.fetch_adcirc_data as fetch_adcirc_data
import gridmap.grid_to_station_maps as grid_to_station_maps

from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities

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
    print(type_pair)
    cast_type = type_pair.split('_')[1]
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

    main_config = utilities.init_logging(subdir=None, config_file='./config/main.yml')

    # Grab the args
    map_file=args.map_file
    url = args.url
    return_sample_min=15
    ensemble=args.ensemble

    stoptime=args.stoptime
    ndays=args.ndays

    if url is None:
        utilities.log.error('For URL-based approaches template url must be specified')
        sys.exit(1)

    if map_file is None: 
        utilities.log.error('map_file must be specified')
        sys.exit(1)

    # Get the IO basic s setep
    utilities.log.info("Product Level Working in {}.".format(os.getcwd()))
    rootdir=io_utilities.construct_base_rootdir(main_config['DEFAULT']['RDIR'], base_dir_extra='')

    if args.sources:
         print('Return list of sources')
         return SOURCES
         sys.exit(0)
    data_source = args.data_source
    data_product = args.data_product

    if data_source.upper() in SOURCES:
        utilities.log.info('Found selected data source {}'.format(data_source))
    else:
        utilities.log.error('Invalid data source {}'.format(data_source))
        sys.exit(1)

    # Setup times range tuple (These could be times OR advisories)
    endtime=args.stoptime
    starttime = genurls.construct_starttime_from_offset(endtime, ndays)

    utilities.log.info('Selected time range is {} to {}, ndays is {}'.format(starttime,endtime,args.ndays))
    time_range=(starttime,endtime)

    #time_start=time_stop+dt.timedelta(days=args.ndays) # How many days BACK
    #starttime=dt.datetime.strftime(time_start, dformat)
    #endtime=dt.datetime.strftime(time_stop, dformat)
    #utilities.log.info('Selected time range is {} to {}, ndays is {}'.format(starttime,endtime,args.ndays))
    #time_range=(starttime,endtime)

    # Construct the URL list for processing: Note this can be used for Hurricane urls as well
    # A little code trickery here to get grid_name
    rpl = genurls.generate_urls_from_times(url=url, timeout=stoptime, ndays=ndays)
    urls = rpl.build_url_list_from_template_url_and_offset(ensemble=args.ensemble)
    grid_name = fetch_adcirc_data.grab_gridname_from_url(urls)
    print(urls)

    # Fetch the grid specific station list. No controls required
    station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=grid_name, mapfile=map_file)
    print(station_file)

    rpl = get_adcirc_stations.get_adcirc_stations(source='ASGS', product=args.data_product,
                station_list_file=station_file)

    # Ensures URLs are desired fort type. 
    urls=get_adcirc_stations.convert_urls_to_61style(urls)
    print(urls)

    # Fetch best resolution and no resampling
    data_adc,meta_adc=rpl.fetch_station_product(urls, return_sample_min=return_sample_min ) # Always grab 15mn sampling

    # Get extra info for amending output filenames
    sitename=rpl.sitename
    gridname=rpl.gridname
    ensemble=rpl.ensemble
    # Lastly fetch the final time of the actual data.
    print(data_adc.index)
    earliest_real_time = min(data_adc.index.tolist()).strftime(dformat)
    latest_real_time = max(data_adc.index.tolist()).strftime(dformat)
    
    # Look inside meta and see what kind of data this is NOWCAST/FORECAST
    cast_type = find_cast_status(meta_adc)

    # Output the data files for subsequent ingesting
    df_adcirc_data = format_data_frames(data_adc) # WIll reformat time index to strings
    adcirc_metadata=sitename.upper()+'_'+ensemble.upper()+'_'+grid_name.upper()+'_'+cast_type+'_'+endtime.replace(' ','T')+'_'+earliest_real_time.replace(' ','T')+'_'+latest_real_time.replace(' ','T')
    try:
        dataf=io_utilities.write_csv(df_adcirc_data, rootdir=rootdir,subdir='',fileroot='adcirc_stationdata',iometadata=adcirc_metadata)
        metaf=io_utilities.write_csv(meta_adc, rootdir=rootdir,subdir='',fileroot='adcirc_stationdata_meta',iometadata=adcirc_metadata)
        utilities.log.info('ADCIRC data has been stored {},{}'.format(dataf,metaf))
    except Exception as e:
        utilities.log.error('Error: ADCIRC: Failed Write {}'.format(e))
        sys.exit(1)

    utilities.log.info('Finished with data source {}'.format(data_source))
    utilities.log.info('Finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--ndays', action='store', dest='ndays', default=-2, type=int,
                        help='Number of look-back days from stoptime: default -2')
    parser.add_argument('--stoptime', action='store', dest='stoptime', default=None, type=str,
                        help='Desired stoptime YYYY-mm-dd HH:MM:SS')
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
    parser.add_argument('--ensemble', action='store',dest='ensemble', default='nowcast',
                        help='str: Select appropriate ensemble Default is nowcast')
    args = parser.parse_args()
    sys.exit(main(args))
