#!/usr/bin/env python

## This class expects to be execited as part of a loop over local yamls
## We expect but do not enforce a single SOURCE
## For each source we process all products


## Recast this into a Class more suitable for lauching by Pegasus

import os,sys
import shutil
import numpy as np
import pandas as pd
import datetime as dt
import harvester.get_observations_stations as get_obs_stations
from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities
import logging

from pathlib import Path

# --- Global fixed options ----------------------------------------------------------

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

# Map AST source names to actual source names
DICT_SOURCES={'USGS_RIVERS':'USGS_RIVERS','USGS':'USGS','NOAA_STATIONS':'NOAAWEB','NDBC_HISTORIC':'NDBC_HISTORIC','NDBC_BUOYS':'NDBC','CONTRAILS_COASTAL':'CONTRAILS','CONTRAILS_RIVERS':'CONTRAILS'}


# --- Main workflow class ----------------------------------------------------------
class FetchObeservations(object):
    stoptime = None
    starttime = None
    ndays = None
    finalDIR = None
    contrails_auth = None
    map_source_file = None
    sampling_min = None
    nwindow = None
    extra_time_stop = 120 # If selecting prediction data try to look forward by 120 hours.
    smooth = None
    data_source = None # In contrast to general hervester, we expect only a single source for these runs
    data_products = None # A user may have multiple products specified, eg wind speed, air pressure, water level etc

    # --- Init ---------------------------------------------------------------------
    def __init__(self, starttime, stoptime, source_config, finalDIR, sampling_min, contrails_auth, smooth, nwindow=11): # , dagfile="workflow.yml"):
        """ 
            starttime: format '%Y-%m-%d %H:%M:%S': Eg --stoptime "2025-08-01 00:00:00".
            stoptime: format '%Y-%m-%d %H:%M:%S': Eg --stoptime "2025-08-01 00:00:00".
            contrails_auth: File containined login Key if wanting to fetch Contrails data else not required.
            smooth: Boolean. Returns data of rolling averaged (centered) windows of width nwindow.
            sampling_min : The time frequency (minutes) for the output timeseries. If req, data are linearly inteprolated but only fills single gaps.
            finalDIR: Location to write all data files.
            source_config: A dict within which are specified the desired list of stations and one or more data products.
        """ 
        self.starttime = starttime
        self.stoptime = stoptime
        self.finalDIR = finalDIR
        self.sampling_min = sampling_min
        self.contrails_auth = contrails_auth
        self.do_smooth = smooth
        self.nwindow = nwindow 
        self.source_config = source_config
        data_for_write = dict()

    def get_timerange(self):
        return (self.starttime,self.stoptime)

    def summary_data(self):
        print(f'Starttime={self.starttime}, Stoptime={self.stoptime}, Number of days value={self.ndays}')

    def construct_root_filenames(self, data_source):
        """
        Define (partial) metadata for the current source
        """
        output_fileroot=f'{data_source.lower()}_stationdata'
        output_metafileroot=f'{data_source.lower()}_stationdata_meta'
        return output_fileroot,output_metafileroot

    def list_sources(self):
        print(f'Current Sources: {list(self.source_config.keys())}')
        return list(self.source_config.keys())

    def list_products(self, insrc):
        print(f'{list(self.source_config[insrc]["SOURCES"].values())}')

    def fetch_data_sources(self): 
        """
        Users may have more than one source in a yaml but will lose parallelsm if they do
        """
        data_for_write = dict()
        print(f'Full dict {self.source_config}')

        for data_source_in in self.list_sources():
            dataf,metaf = self.fetch_data_products(data_source_in) #  self.source_config[data_source_in]['SOURCES'].values())
            data_for_write[data_source_in]=[dataf,metaf]
        return data_for_write

    def fetch_data_products(self, data_source_in):
        """
            For the input SOURCEs, grab data for all of the requested products. For increased parallelism, 
            the user should only specify a singe SOURCE and a single PRODUCT per input yaml.

            If smoothing we still need the raw data but fetch it using maximum av ailable resolution 
            via sampling_min 0. We also must must remove daa with lots of missingness - also update meta
        """

        raw_sampling_min = self.sampling_min if not self.do_smooth else 0
        outdir=self.finalDIR

        data_source=DICT_SOURCES[data_source_in]

        # This one needs the NOAA_STATION word
        source_products = list(self.source_config[data_source_in]['SOURCES'].values())

        station_file=self.source_config[data_source_in]['STATION_FILE']
        utilities.log.debug(f' station file {station_file}')
        data_source_short = 'CONTRAILS' if 'CONTRAILS' in data_source_in else data_source
        output_fileroot,output_metafileroot = self.construct_root_filenames(data_source_short)
        utilities.log.debug(f'File rootnames are {output_fileroot} and {output_metafileroot}')

        endt=self.stoptime.replace(' ','T') # For filename metadata

        datafs = list()
        metafs=list()
        for data_product in source_products:
            print(f'Preparing to fetch {data_product} from {data_source_in}')

            metadata = self.construct_metadata(data_product, endt)
            obs_get = get_obs_stations.get_obs_stations(source=data_source_short.upper(), product=data_product,
                contrails_yamlname=self.contrails_auth, knockout_dict=None, station_list_file=station_file)

            ## Expose option to SMOOTH data to the desired window width
            ## If smooth is selected choose the full resolution to get raw data, smooth and then reset freq
            try: # Do this in case we got no data from any station
                df_data,meta=obs_get.fetch_station_product(self.get_timerange(), return_sample_min=raw_sampling_min)
                df_data.index = pd.to_datetime(df_data.index) # Do this now and it propogates to the subsequent daasets
                if self.do_smooth:
                    utilities.log.info('Smoothed data requested with sampling_min {sampling_min}')
                    data_thresholded = obs_get.remove_missingness_stations(data, max_nan_percentage_cutoff=100) # (maxmum percentage of allowable nans)
                    meta_list = set(data_thresholded.columns.tolist()).intersection(meta.index.to_list())
                    meta = meta.loc[meta_list]
                    data_smoothed = obs_get.fetch_smoothed_station_product(data_thresholded, return_sample_min=sampling_min, window=nwindow)
                    df_data = data_smoothed # format_data_frames(data_smoothed, data_product) 
                    metadata=f'{metadata}_smoothed'
                dataf,metaf = self.write_data_todisk(outdir,data_source_short,df_data,meta,metadata,output_fileroot,output_metafileroot)
                datafs.append(dataf)
                metafs.append(metaf)
            except Exception as e: # IndexError as e:
                print(f'Failed data products {data_product} for source {data_source_in}. Skip e={e}')
                pass
            utilities.log.info(f'Finished with data source {data_source_in}')
            return dataf, metaf 

    def construct_metadata(self, data_product, endt):
        """   
        Construct remaining file metadata and include the stoptime values needed for DB uploads
        Complexity caused by this provider segregating rivers from coastal data
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

    def write_data_todisk(self, outdir,data_source,data,meta,metadata,output_fileroot,output_metafileroot): 
        try:
            dataf=io_utilities.write_csv(data, rootdir=outdir,subdir='',fileroot=output_fileroot,iometadata=metadata)
            print(f'{data_source.upper()} data has been stored {dataf}')
            metaf=io_utilities.write_csv(meta, rootdir=outdir,subdir='',fileroot=output_metafileroot,iometadata=metadata)
            print(f'{data_source.upper()} metadata has been stored {metaf}')
        except Exception as e:
            print(f'Error:{data_source} Failed Write {e}')
            sys.exit(1)
        return dataf,metaf

# --- Run an example from the CLI --------------------------------------------------

def main(args):
    """
    Generally we anticipate inputting a STOPTIME
    Then the STARTTIME is ndays on the past
    """

    main_config = utilities.init_logging(subdir=None,config_file=os.path.join(os.path.dirname(__file__),'./config','main.yml'))

    print("Product Level Working in {}.".format(os.getcwd()))

    BASE_DIR = Path(".").resolve()
    print(f' BASE_DIR {BASE_DIR}')
    
    # Maybe mkdir ./outputs ?

    utilities.log.info("Product Level Working in {}.".format(os.getcwd()))

## Select where to write these files

    finalDIR=io_utilities.construct_base_rootdir(args.finalDIR, base_dir_extra=None)

    ## Return frequency of the data
    sampling_min=args.sampling_min

    ## If choosing smooth only smoothed data are returned
    do_smooth = args.smooth
    nwindow = args.nwindow if do_smooth else None

    print(f'Requested sampling frequency {sampling_min}')
    print(f'Output directory specified to be {finalDIR}')
    print(f'Output user directive for smoothing {do_smooth}')

    contrails_auth=args.contrails_auth
    sampling_min = args.sampling_min
    map_source_file = args.map_source_file
    #finalDIR = args.finalDIR

    # Get local list of sources/products to process
    source_config = utilities.load_config(args.map_source_file)['SOURCEMAP']

    # Get stoptime and infer starttime
    time_stop = dt.datetime.now() if args.stoptime is None else dt.datetime.strptime(args.stoptime,dformat)
    ndays = args.ndays
    time_start=time_stop+dt.timedelta(days=args.ndays) # How many days BACK
    starttime=dt.datetime.strftime(time_start, dformat)
    endtime=dt.datetime.strftime(time_stop, dformat)

    ## Setup Job 

    obs = FetchObeservations(starttime, endtime, source_config, finalDIR, sampling_min, contrails_auth, do_smooth, nwindow) 
    obs.summary_data()

    time_range=obs.get_timerange() # (starttime,endtime)

    # Ths is out of place metadata are used to augment filename
    endt=endtime.replace(' ','T')

    # Process the sources/products from the map_source_file
    data_sources = obs.list_sources() # list(source_config.keys())

    # Fetch the data
    data_for_write = obs.fetch_data_sources() # Aggregate of the data and where they should be written
    print(data_for_write)

    ## Now write out the data 
#   def write_data_todisk(self, outdir,data_source,data,meta,metadata,output_fileroot,output_metafileroot):

    print('Finished') 

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
    parser.add_argument('--sampling_min', action='store', dest='sampling_min', default=15, type=int,
                        help='Returned data sampling frequency: default 15 (min)')
    parser.add_argument('--smooth', action='store_true', default=False,
                        help='Boolean: Will inform Harvester to return smoothed data')
    parser.add_argument('--nwindow', action='store', dest='nwindow', default=11, type=int,
                        help='Rolling window width in steps of sapling_min units. If smooth is selecte: default 11')

    args = parser.parse_args()
    sys.exit(main(args))
