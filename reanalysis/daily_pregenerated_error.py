#!/usr/bin/env python
#
# Methods to read the input time series (time x stations) for a single year
# For that year, compute a lowpass FFT.
# Resample the data at each day (00Z) and construct a new, daily, stationSummaryAves.csv file 
#    for building daily interpolation field corrections.
#    This choice of resampling may cause the first day of the year to be removed as it usually starts at 01Z.
#
# Lastly, since it has been useful in the past, construct a set of means for comparison
#

# This code does very little but is retained to keep our reanalysis pipeline intact
# Basically, we have already estimated all the 43 years of errors for each station. So all we need to do here
# is write the files into the appropaiate errorfield subdirectory. Also, write out all 4 data sets for the desired timemark:
# SLR,MID,VLF,SUM

import os,sys
import numpy as np
import pandas as pd
import time as tm
import seaborn as sns
from pandas.tseries.frequencies import to_offset
import json
from scipy import signal
from scipy.signal import butter, lfilter, savgol_filter
import datetime as dt

import plot_daily_lowpass_sampled_error as plot_daily
from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities

def generate_errorfiles(iometa, df_all, df_meta, outputdir=None, fileroot='stationSummaryAves'):
    """
    """
    subdir='errorfield'
    datadict = dict()
    print(f'Total errorrows {df_all}')
    for index, df in df_all.iterrows():
        print(f' ITERATE {index}')
        print(f' df {df}')
        utilities.log.info(f' Index data {index}')
        utilities.log.info(f' IOMETA {iometa}')
        utilities.log.info(f' IOMETA keys {iometa.keys()}')
        utilities.log.info(f' IOMETA[index] {iometa[index]}')
        df = df.to_frame()
        utilities.log.info(f' df {df}')
        df.index.name='STATION'
        df.columns=['VAL'] 
        metadata=iometa[index]
        ##df.index = df.index.astype('int64')    
        df_merged=df_meta.join(df)
        df_merged=df_merged[['LAT','LON','NAME','VAL']]
        utilities.log.info(f'Data shape after nan drop: {df_merged.shape}')
        outfilename=io_utilities.write_csv(df_merged,rootdir=outputdir,subdir=subdir,fileroot=fileroot,iometadata=metadata)
        datadict[iometa[index]]=outfilename
        df_merged.to_csv(outfilename)
        print(f'Wrote out the summary files {metadata} and file {outfilename}')
    return datadict

# The means are proper means for the month./week. The offset simply changes the index underwhich it will be stored
# The daily means index at starting week +3 days.
def station_level_means(df_obs, df_adc, df_err, station):
    dfs = pd.DataFrame()
    dfs['OBS']=df_obs[station]
    dfs['ADC']=df_adc[station]
    dfs['ERR']=df_err[station]
    #
    ###### Doesnt work when spanning multiple yearsdfs['Year'] = dfs.index.year
    dfs['Month'] = dfs.index.month
    dfs['Hour'] = dfs.index.hour
    #
    data_columns = ['OBS','ADC','ERR']
    # Resample to monthly frequency, aggregating with mean
    dfs_monthly_mean = dfs[data_columns].resample('MS').mean() #MS restarts meaning from the start not the end
    dfs_monthly_mean.index = dfs_monthly_mean.index + to_offset("15d")
    dfs_weekly_mean = dfs[data_columns].resample('W').mean()
    # Get rolling
    dfs_7d = dfs[data_columns].rolling(7, center=True).mean()
    return dfs, dfs_weekly_mean, dfs_monthly_mean, dfs_7d

def dict_to_data_frame(dataDict, src):
    stations = list(dataDict.keys())
    dt = [dataDict[x][src]['TIME'] for x in dataDict.keys()]
    #dx = [dataDict[x][src]['WL'] for x in dataDict.keys()]
    indexes = dt[0]
    df=pd.DataFrame(indexes)
    df.columns=['TIME']
    for station in stations:
        df[station]=dataDict[station][src]['WL']
    df.set_index('TIME',inplace=True)
    df.index = pd.to_datetime(df.index)
    return df

# OrderedDicts ?
def fetch_data_metadata(f, meta):
    with open(meta, 'r') as fp:
        try:
            metaDict = json.load(fp)
        except OSError:
            utilities.log.error("Could not open/read file {}".format(meta)) 
            sys.exit()
    # Timeseries data
    with open(f, 'r') as fp1:
        try:
            dataDict = json.load(fp1)
        except OSError:
            utilities.log.error("Could not open/read file {}".format(meta))
            sys.exit()
    return dataDict, metaDict

def main(args):
    t0 = tm.time()

    #assumes the existance of a proper main.yml to get IO information
    if args.main_yamlname is None:
        main_yamlname=os.path.join(os.path.dirname(__file__), './config', 'main.yml')
    else:
        main_yamlname=args.main_yamlname
    config = utilities.init_logging(subdir=None, config_file=main_yamlname, log_file_metadata='daily')
    utilities.log.info('Selected main_yamlfile of {}'.format(main_yamlname))

    if not args.input_directory:
        utilities.log.error('Need input_directory on command line: --input_directory <input_directory>')
        return 1
    topdir = args.input_directory.strip()

    if not args.output_directory:
        utilities.log.error('Need output_directory on command line: --input_directory <input_directory>')
        return 1

    outputdir=io_utilities.construct_base_rootdir(args.output_directory.strip(), base_dir_extra='')

    # Define a generic time range for any year: Maybe want to include the flank
    inyear = args.inyear
    timein = '-'.join([inyear,'01','01'])
    timeout = '-'.join([inyear,'12','31','18'])

    # Grab filenames for generation of the error data
    yearly_file_names = io_utilities.read_json_file(f'{topdir}/{args.yearly_file_properties}')
    f = yearly_file_names['adc_obs_error_merged'] 
    meta = yearly_file_names['obs_wl_metadata_json']

    dataDict, metaDict = fetch_data_metadata(f, meta)

    stations = list(dataDict.keys()) # For subsequent naming - are we sure order is maintained?
    utilities.log.info(f'Input data chosen range is {timein}, {timeout}')
    utilities.log.info(f'Num of input stations is {len(stations)}')

    # Metadata
    df_meta=pd.DataFrame(metaDict)
    df_meta.index.name='STATION'

    utilities.log.info('Selecting yearly data between {} and {}, inclusive'.format(timein, timeout))
    # Time series data. This ONLY works on compError generated jsons
    df_obs_all = dict_to_data_frame(dataDict, 'OBS').loc[timein:timeout] # Get whole year inclusive
    df_adc_all = dict_to_data_frame(dataDict, 'ADC').loc[timein:timeout]
    df_err_all = dict_to_data_frame(dataDict, 'ERR').loc[timein:timeout]

    ####start = df_err_all.index.min()

    # Fetch the station errror estimations instead of building new filtered data sets
    #
    estimator_dir='/projects/sequence_analysis/vol1/prediction_work/NOAA.Reanalysis/STATION_ML_ESTIMATORS_DETRENDED_OBS-V2-43Years'
    slr_filename=f'{estimator_dir}/SLR/ALL_stations_SLR_estimations.pkl'
    slr_metafilename=f'{estimator_dir}/SLR/ALL_stations_SLR_metadata.pkl'

    mid_filename=f'{estimator_dir}/MIDRANGE_FREQ/ALL_station_bandpass_data.pkl'
    vlf_filename=f'{estimator_dir}/VERYLOWFREQ/ALL_stations_XGBoost_hourly_estimations.pkl'

    # Need to read elsewhere for the metadata 

    hsofs_yearlydir='/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/reanalysis/hsofs.V2/YEARLY-2020'
    mid_metafilename=f'{hsofs_yearlydir}/obs_wl_metadata_hourly_height.pkl'
    vlf_metafilename=f'{hsofs_yearlydir}/obs_wl_metadata_hourly_height.pkl'

    #
    # Read the estimated error models
    #
    df_all_slr = pd.read_pickle(slr_filename)
    df_all_meta_slr = pd.read_pickle(slr_metafilename)

    df_all_mid = pd.read_pickle(mid_filename)
    df_all_meta_mid = pd.read_pickle(mid_metafilename)

    df_all_vlf = pd.read_pickle(vlf_filename)
    df_all_meta_vlf = pd.read_pickle(vlf_metafilename)

    df_all_slr = df_all_slr.loc[timein:timeout].copy()
    df_all_mid = df_all_mid.loc[timein:timeout].copy()
    df_all_vlf = df_all_vlf.loc[timein:timeout].copy()

    ##################################################################################

    # NOTE:, since the SLR data has about 80 stations and the MID/VLF hsofs grid is about 52, 
    # the filename data sets will be of different length. This will get resolved by the subsequent
    # interpolation procedures


    #intersectedStations=set(df_err_all.columns.to_list()).intersection(stations) # Compares data to metadata lists
    #utilities.log.info('Number of intersected stations is {}'.format(len(intersectedStations)))
    #utilities.log.debug('Intersected stations {}'.format(intersectedStations))

    ##
    ## We want to fetch data on 6 hourly intervales instead of days
    ## Per the analyses of Taylor A. and Brian B. 


    # Now pull out daily data 
    utilities.log.info('MANUAL setting of date range')
    #stime=''.join(['1979','-01-01 00:00:00'])
    #etime=''.join(['1979','-05-01 00:00:00'])

    stime=timein
    etime=timeout

    # Switch to doing every 6 hourly moving forward. Do not force starting at 00Z since 1979 starts at 01Z
    starttime = dt.datetime.strptime(stime,'%Y-%m-%d') 
    endtime = dt.datetime.strptime(etime,'%Y-%m-%d-%H') 
    #numDays = (endtime-starttime).days + 1

    utilities.log.info('Specified time ranges are {} and {}'.format(starttime, endtime))

    ##startday=pd.date_range(starttime, periods=numDays) #.values()
    #startday=pd.date_range(starttime, periods=numDays,freq='6H')
    startday=pd.date_range(starttime, endtime,freq='6H')

    ##julianMetadata = startday.strftime('%y-%j').to_list()
    julianMetadata = startday.strftime('%y-%j-%H').to_list()
    utilities.log.info('startdays {}'.format(startday))
    
    # Build metadata for the files
    # Julian dat 
    iometa = dict()
    for index,day,date in zip(range(len(startday)),julianMetadata,startday):
        #iometa[index]='_'.join([day,date.strftime('%Y%m%d%H')])
        iometa[date]='_'.join([day,date.strftime('%Y%m%d%H')])

    # Now process the Rows and build a new datafile for each
    # df_meta and df report stationids as diff types. Yuk.
    # Store the list of filenames into a dict for interpolation processing
    print(f' IOMETA {iometa}')

    intersect = [value for value in startday if value in df_all_slr.index]
    df_all_slr_subselect=df_all_slr.loc[intersect]

    intersect = [value for value in startday if value in df_all_mid.index]
    df_all_mid_subselect=df_all_mid.loc[intersect]

    intersect = [value for value in startday if value in df_all_vlf.index]
    df_all_vlf_subselect=df_all_vlf.loc[intersect]

    datadict_all=dict()
    datadict_all['SLR']=generate_errorfiles(iometa, df_all_slr_subselect, df_all_meta_slr, outputdir=outputdir, fileroot='stationSummarySLR')
    datadict_all['MID']=generate_errorfiles(iometa, df_all_mid_subselect, df_all_meta_mid, outputdir=outputdir,fileroot='stationSummaryMID')
    datadict_all['VLF']=generate_errorfiles(iometa, df_all_vlf_subselect, df_all_meta_vlf, outputdir=outputdir,fileroot='stationSummaryVLF')

    print('Write out runProps.json file')
    outfilesjson = io_utilities.write_dict_to_json(datadict_all, rootdir=topdir,subdir='',fileroot='daily_file_properties',iometadata='')

## Maybe just launch the slurm jobs from here?

    utilities.log.info('Wrote pipeline Dict data to {}'.format(outfilesjson))
    print('Finished generating daily lowpass data files')

if __name__ == '__main__':
    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()
    parser.add_argument('--main_yamlname', action='store',dest='main_yamlname', default=None,
                        help='str: Select appropriate main_yamlname')
    parser.add_argument('--yearly_file_properties', action='store', dest='yearly_file_properties', default=None,
                        help='Dict that maps compute_error keys to filenames for output data sets')
    parser.add_argument('--input_directory', action='store', dest='input_directory', default=None,
                        help='directory for yearly data')
    parser.add_argument('--output_directory', action='store', dest='output_directory', default=None,
                        help='directory for yearly data')
    parser.add_argument('--inyear', action='store', dest='inyear', default=None, help='potentially contains flanks about the main year')
    args = parser.parse_args()
    sys.exit(main(args))
