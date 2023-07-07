#!/usr/bin/env python
#
# A small code to facilitate the gathering of ADCIRC water data using the AST:harvester classes
# The expectation is this code is executed by a bash script within a cron 
#
# This wrapper only fetches ADCIRC data using the URL Template approach
# Only a fort61_style update is performed for this data execution
#
# This approach expects to fetch urlsa based on the url_framework.yml

#python get_adcirc_stations.py --url "http://tds.renci.org/thredds/dodsC/2021/al09/11/hsofs/hatteras.renci.org/hsofs-al09-bob/nhcOfcl/fort.61.nc" --data_source 'TDS' --timeout 8  --fort63_style



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
    data_list = df['CAST'].tolist()
    cast_type = list(set(data_list))[0] # fetch_stations has already ensured these will all be the same grid and type
    utilities.log.info(cast_type)
    return cast_type.upper()

dformat='%Y-%m-%d %H:%M:%S'
GLOBAL_TIMEZONE='gmt' # Every source is set or presumed to return times in the zone
PRODUCT='water_level'

# Currently supported sources
SOURCES = ['TDS']

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
    config_name=args.config_file
    return_sample_min=15
    instance_name=args.instance_name
    grid_name=args.grid_name
    ensemble=args.ensemble

    stoptime=args.stoptime
    ndays=args.ndays

    # Somewhat unusual extras
    hur_source=args.hurricane_source
    hur_year=args.hurricane_year

    if instance_name is None or ensemble is None or grid_name is None:
        utilities.log.error('For YAML-based approaches all of instance,grid_name,ensemble must be specified')
        sys.exit(1)

    if map_file is None or config_name is None:
        utilities.log.error(' Both the map_file and config_file must be specified')
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

    # Setup times/advisories range tuple

    endtime=args.stoptime
    starttime = genurls.construct_starttime_from_offset(endtime, ndays) 

    utilities.log.info('Selected time range is {} to {}, ndays is {}'.format(starttime,endtime,args.ndays))
    time_range=(starttime,endtime)

    # Fetch the grid specific station list. No controls required
    station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.grid_name, mapfile=map_file)

    # Construct the URL list for processing: Note this can be used for Hurricane urls as well
    rpl = genurls.generate_urls_from_times(timeout=stoptime, ndays=ndays, grid_name=grid_name, instance_name=instance_name, config_name=config_name, hurricane_yaml_year=hur_year,hurricane_yaml_source=hur_source)
    urls = rpl.build_url_list_from_yaml_and_offset(ensemble=args.ensemble)
    print(urls)
    print(station_file)

    rpl = get_adcirc_stations.get_adcirc_stations(source='TDS', product=args.data_product,
                station_list_file=station_file)

    timemark=rpl.timemark

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
    df_adcirc_data = format_data_frames(data_adc)

    try:
        tone=dt.datetime.strptime(timemark,'%Y%m%d%H')
        timemark=dt.datetime.strftime(tone, dformat)
        timemark=timemark.replace(' ','T')
    except Exception as e:
        pass

    adcirc_metadata=sitename.upper()+'_'+ensemble.upper()+'_'+grid_name.upper()+'_'+cast_type+'_'+timemark+'_'+earliest_real_time.replace(' ','T')+'_'+latest_real_time.replace(' ','T')

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
                        help='Number of look-back days from stoptime : default -2')
    parser.add_argument('--stoptime', action='store', dest='stoptime', default=None, type=str,
                        help='Desired stoptime YYYY-mm-dd HH:MM:SS')
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--data_source', action='store', dest='data_source', default='TDS', type=str,
                        help='choose supported data source: default = TDS')
    parser.add_argument('--data_product', action='store', dest='data_product', default='water_level', type=str,
                        help='choose supported data product: default is water_level')
    parser.add_argument('--map_file', action='store', dest='map_file', default=None, type=str,
                        help='Location of the grid_to_stationfile_maps.yml data')
    parser.add_argument('--config_file', action='store', dest='config_file', default=None, type=str,
                        help='Location of the url_framework.yml data')
    parser.add_argument('--instance_name', action='store', dest='instance_name', default=None,
                        help='String: instance name')
    parser.add_argument('--hurricane_source', action='store',dest='hurricane_source', default=None,
                        help='str: Only needed for Hurricanes AND if using YAML to specify urls')
    parser.add_argument('--hurricane_year', action='store',dest='hurricane_year', default=None,
                        help='str: Only needed for Hurricanes AND if using YAML to specify urls')
    parser.add_argument('--grid_name', action='store',dest='grid_name', default=None,
                        help='str: Select appropriate grid_name')
    parser.add_argument('--ensemble', action='store',dest='ensemble', default='nowcast',
                        help='str: Select appropriate ensemble Default is nowcast')

    args = parser.parse_args()
    sys.exit(main(args))
