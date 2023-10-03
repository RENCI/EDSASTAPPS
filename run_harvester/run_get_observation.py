#!/usr/bin/env python
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
    prod_name = PRODUCT[product].upper() if product in PRODUCT.keys() else 'NONAME'
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
        'river_flow_volume': 'river_flow_volume',
        'river_stream_elevation':'water_level',
        'wind_speed':'wind_speed'
         }

# Currently supported sources
SOURCES = ['NOAAWEB', 'NOAA','CONTRAILS','NDBC']

##
## Run stations
##

def main(args):
    """
    Generally we anticipate inputting a STOPTIME
    Then the STARTTIME is ndays on the past
    """

    main_config = utilities.init_logging(subdir=None,config_file=os.path.join(os.path.dirname(__file__),'./config','main.yml'))

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

    # Setup times range tuple
    if args.stoptime is not None:
        time_stop=dt.datetime.strptime(args.stoptime,dformat)
    else:
        time_stop=dt.datetime.now()
    time_start=time_stop+dt.timedelta(days=args.ndays) # How many days BACK
    starttime=dt.datetime.strftime(time_start, dformat)
    endtime=dt.datetime.strftime(time_stop, dformat)
    utilities.log.info('Selected time range is {} to {}, ndays is {}'.format(starttime,endtime,args.ndays))

    time_range=(starttime,endtime)

    # metadata are used to augment filename
    # Add data_product
    endt=endtime.replace(' ','T')
    metadata=f'{args.data_product}_{endt}' 

    in_station_list = args.station_list
    if in_station_list is not None:
        utilities.log.info('User supplied station list found {}'.format(in_station_list))
    if data_source.upper()=='NOAA':
        utilities.log.info('Preparing for a NOAA fetch')
        time_range=(starttime,endtime) # Can be directly used by NOAA 
        station_file=os.path.join(os.path.dirname(__file__),'./supporting_data','noaa_stations.csv') if in_station_list is None else in_station_list
        output_fileroot='noaa_stationdata'
        output_metafileroot='noaa_stationdata_meta'
    if data_source.upper()=='NOAAWEB':
        utilities.log.info('Preparing for a NOAAWEB fetch')
        time_range=(starttime,endtime) # Can be directly used by NOAA 
        station_file=os.path.join(os.path.dirname(__file__),'./supporting_data','noaa_stations.csv') if in_station_list is None else in_station_list
        output_fileroot='noaa_stationdata'
        output_metafileroot='noaa_stationdata_meta'

    elif data_source.upper()=='CONTRAILS':
        utilities.log.info('Preparing for a CONTRAILS fetch')
        contrails_config = args.config_name # utilities.load_config(os.path.join(os.path.dirname(__file__),'./secrets','contrails.yml'))['DEFAULT']
        utilities.log.info('Got Contrails access information {}'.format(contrails_config))
        if data_product=='river_water_level' or data_product=='river_stream_elevation':
            station_file=os.path.join(os.path.dirname(__file__),'./supporting_data','contrails_stations_rivers.csv') if in_station_list is None else in_station_list
            meta='RIVERS'
        else:
            station_file=os.path.join(os.path.dirname(__file__),'./supporting_data','contrails_stations_coastal.csv') if in_station_list is None else in_station_list
            meta='COASTAL'
        metadata=f'{metadata}_{meta}'
        output_fileroot='contrails_stationdata'
        output_metafileroot='contrails_stationdata_meta'
    elif data_source.upper()=='NDBC':
        utilities.log.info('Preparing for a NDBC fetch')
        time_range=(starttime,endtime) # Can be directly used by NDBC
        station_file=os.path.join(os.path.dirname(__file__),'./supporting_data','ndbc_buoys.csv') if in_station_list is None else in_station_list
        output_fileroot='ndbc_stationdata'
        output_metafileroot='ndbc_stationdata_meta'
    else:
        utilities.log.error('Failed: Only NDBC, NDBC_HISTORIC, NOAAWEB, NOAA or CONTRAILS currently supported')
        sys.exit(1)

    obs = get_obs_stations.get_obs_stations(source=data_source.upper(), product=args.data_product,
            contrails_yamlname=args.config_name, knockout_dict=None, station_list_file=station_file)

    # Get data at highest resolution. Return at 15min intervals
    data,meta=obs.fetch_station_product(time_range, return_sample_min=15, interval='None' )

    df_data = format_data_frames(data, data_product) # Melt data to support later database updates

    # Output
    try:
        dataf=io_utilities.write_csv(df_data, rootdir=rootdir,subdir='',fileroot=output_fileroot,iometadata=metadata)
        metaf=io_utilities.write_csv(meta, rootdir=rootdir,subdir='',fileroot=output_metafileroot,iometadata=metadata)
        utilities.log.info('OBS data has been stored {},{}'.format(dataf,metaf))
    except Exception as e:
        utilities.log.error('Error:{} Failed Write {}'.format(args.data_source,e))
        sys.exit(1)

    utilities.log.info('Finished with data source {}'.format(data_source))

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--ndays', action='store', dest='ndays', default=-2, type=int,
                        help='Number of look-back days from stoptime (or now): default -2')
    parser.add_argument('--stoptime', action='store', dest='stoptime', default=None, type=str,
                        help='Desired stoptime YYYY-mm-dd HH:MM:SS. Default=now')
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--data_source', action='store', dest='data_source', default=None, type=str,
                        help='choose supported data source (case independant) eg NOAA, NOAAWEB or CONTRAILS')
    parser.add_argument('--data_product', action='store', dest='data_product', default=None, type=str,
                        help='choose supported data product eg river_water_level: Only required for Contrails')
    parser.add_argument('--station_list', action='store', dest='station_list', default=None, type=str,
                        help='Choose a non-default location/filename for a stationlist')
    parser.add_argument('--config_name', action='store', dest='config_name', default=None, type=str,
                        help='Choose a non-default contrails auth config_name')
    parser.add_argument('--finalDIR', action='store', dest='finalDIR', default=None,
                        help='String: Custom location for the output dicts, PNGs and logs')
    args = parser.parse_args()
    sys.exit(main(args))
