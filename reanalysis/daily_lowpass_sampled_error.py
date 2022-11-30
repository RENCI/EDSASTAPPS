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
# This code is a rewrite of earlier reanalysis code. Mostly changed to use the new AST data and formats. 
#

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

##
## Define different choices for bandpass filtering
##

def butter_lowpass_filter(df_data, filterOrder=10, numHours=200):
    """
    Butterworth lowpass filter.
    Note. Need to drop nans from the data set
    """
    numHoursCutoff=numHours
    cutoff = 2/numHoursCutoff #Includes nyquist adjustment
    #ddata=df[station]
    ddata = df_data
    print('cutoff={} Equiv number of hours is {}'.format(cutoff, 2/cutoff))
    print('filterOrder is {}'.format(filterOrder))
    sos = signal.butter(filterOrder, cutoff, 'lp', output='sos',analog=False)
    datalow = signal.sosfiltfilt(sos, ddata)
    #print('ddata {}'.format(ddata))
    #print('LP data {}'.format(datalow))
    return datalow

def fft_lowpass(signal, lowhrs):
    """
    Performs a low pass filer on the series.
    low and high specifies the boundary of the filter.

    NOTE for high freq we always take the highest possible
    which is just the sampling freq
    """
    low = 1/lowhrs
    high = 1/signal.shape[0]
    print('FFT: low {}, high {}'.format(low,high))
    if len(signal) % 2:
        result = np.fft.rfft(signal, len(signal))
    else:
        result = np.fft.rfft(signal)
    freq = np.fft.fftfreq(len(signal))[:len(signal) // 2 + 1]
    factor = np.ones_like(freq)
    factor[freq > low] = 0.0
    sl = np.logical_and(high < freq, freq < low)
    a = factor[sl]
    # Create float array of required length and reverse.
    a = np.arange(len(a) + 2).astype(float)[::-1]
    # Ramp from 1 to 0 exclusive.
    a = (a / a[0])[1:-1]
    # Insert ramp into factor.
    factor[sl] = a
    result = result * factor
    ###print('Freqs {}'.format(result))
    ###print('FREQ TYPEW {}'.format(type(result)))
    return np.fft.irfft(result, len(signal))

# The means are proper means for the month./week. The offset simply changes the index underwhich it will be stored
# The daily means index at starting week +3 days.
def station_level_means(df_obs, df_adc, df_err, station):
    dfs = pd.DataFrame()
    dfs['OBS']=df_obs[station]
    dfs['ADC']=df_adc[station]
    dfs['ERR']=df_err[station]
    #
    ###### Doesn;t work when spanning multiple yearsdfs['Year'] = dfs.index.year
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
    timeout = '-'.join([inyear,'12','31'])

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

    start = df_err_all.index.min()

    # FFT Lowpass each station for entire range time. Then, extract values for all stations every day at 00Z
    upshift=0
    hourly_cutoff= 24 # 168 # 48 # 6 # 168 #48
    cutoff = hourly_cutoff+upshift
    utilities.log.info('FFT hourly_cutoff {}, actual_cutoff {}'.format(hourly_cutoff,cutoff))

    intersectedStations=set(df_err_all.columns.to_list()).intersection(stations) # Compares data to metadata lists
    utilities.log.info('Number of intersected stations is {}'.format(len(intersectedStations)))
    #utilities.log.debug('Intersected stations {}'.format(intersectedStations))

    # Perform FFT for each station over the entire time range
    df_err_all_lowpass=pd.DataFrame(index=df_err_all.index)

    fftAllstations=dict()
    io_utilities.write_pickle(df_err_all, rootdir=outputdir,subdir='',fileroot='df_err_all',iometadata='')

    # Remove stations from intersectedStation that may not complete the FFT step
    lostStations=list()
    for station in intersectedStations:
        try:
            print('Process station {}'.format(station))
            stationName = df_meta.loc[station]['NAME']
            df_fft=pd.DataFrame()
            df_low = pd.DataFrame()
            utilities.log.info('Process cutoff {} for station {}'.format(cutoff,station))
            df_temp = df_err_all[station].dropna()
            df_low['low'] = fft_lowpass(df_temp,lowhrs=cutoff)
            df_low.index=df_temp.index 
            df_fft[str(cutoff)]=df_low['low']
            df_err_all_lowpass[station]=df_fft[str(cutoff)]
        except:
            utilities.log.info('FFT failed for station {}'.format(station))
            utilities.log.info('Removing a failed station from the FFT data set {}'.format(station))
            lostStations.append(station)

    utilities.log.info(f'Lost stations {lostStations}')
    #intersectedStations = [x for x in intersectedStations if x not in lostStations] 
    io_utilities.write_pickle(df_err_all_lowpass, rootdir=outputdir,subdir='',fileroot='df_err_all_lowpass',iometadata='')

    # Now pull out daily data 
    utilities.log.info('MANUAL setting of date range')
    #stime=''.join(['1979','-01-01 00:00:00'])
    #etime=''.join(['1979','-05-01 00:00:00'])

    stime=timein
    etime=timeout

    starttime = dt.datetime.strptime(stime,'%Y-%m-%d') 
    endtime = dt.datetime.strptime(etime,'%Y-%m-%d') 
    numDays = (endtime-starttime).days + 1

    utilities.log.info('Specified time ranges are {} and {}'.format(starttime, endtime))

    startday=pd.date_range(starttime, periods=numDays) #.values()
    julianMetadata = startday.strftime('%y-%j').to_list()
    utilities.log.info('startdays {}'.format(startday))
    
    # Build metadata for the files
    # Julian dat 
    iometa = dict()
    for index,day,date in zip(range(len(startday)),julianMetadata,startday):
        #iometa[index]='_'.join([day,date.strftime('%Y%m%d%H')])
        iometa[date]='_'.join([day,date.strftime('%Y%m%d%H')])

    # startday is in yyyy-mm-dd lowpass includes times

    utilities.log.info('df_err_all times {}'.format(df_err_all.index))
    utilities.log.info('df_err_lowpass times {}'.format(df_err_all_lowpass.index))
    intersect = [value for value in startday if value in df_err_all_lowpass.index] 

    utilities.log.info('Residual data: intersect list {}'.format(len(intersect)))
    df_err_all_lowpass_subselect=df_err_all_lowpass.loc[intersect]

    # Now process the Rows and build a new datafile for each
    # df_meta and df report stationids as diff types. Yuk.
    # Store the list of filenames into a dict for interpolation processing

    subdir='errorfield'
    datadict = dict()
    for index, df in df_err_all_lowpass_subselect.iterrows():
        utilities.log.info(f' Index {index}')
        df = df.to_frame()
        utilities.log.info(f' df {df}')
        df.index.name='STATION'
        df.columns=['fft'] 
        metadata=iometa[index]
        ##df.index = df.index.astype('int64')    
        df_merged=df_meta.join(df)
        ##df_merged.drop('NAME',axis=1, inplace=True)
        df_merged=df_merged[['LAT','LON','NAME','STATE','fft']]
        # utilities.log.info(f'Full merged data set {df_merged}')
        # We do not need the STATE for the subsequent analysis
        # df_merged[['LAT','LON','NAME','fft']].dropna(inplace=True) # Cannot pass Nans to the interpolation system
        utilities.log.info(f'Data shape after nan drop: {df_merged.shape}')
        outfilename=io_utilities.write_csv(df_merged,rootdir=outputdir,subdir=subdir,fileroot='stationSummaryAves',iometadata=metadata)
        datadict[iometa[index]]=outfilename
        df_merged.to_csv(outfilename)

    print('Write out runProps.json file')
    outfilesjson = io_utilities.write_dict_to_json(datadict, rootdir=topdir,subdir='',fileroot='daily_file_properties',iometadata='')

## Maybe just launch the slurm jobs from here?

# Run the plot pipeline ASSUMES outputdir has been created already
    sns.set(rc={'figure.figsize':(11, 4)}) # Set gray background and white gird
    for station in stations:
        dfs, dfs_weekly_mean, dfs_monthly_mean, dfs_7d = station_level_means(df_obs_all, df_adc_all, df_err_all, station)
        start=dfs.index.min().strftime('%Y-%m')
        end=dfs.index.max().strftime('%Y-%m')
        stationName = df_meta.loc[station]['NAME']
        plot_daily.make_plot(start, end, station, 'ERR', stationName, dfs, dfs_7d, dfs_weekly_mean, dfs_monthly_mean, outputdir) 
        #make_lowpass_plot no longer works as is
        #plot_daily.make_lowpass_plot(start, end, lowpassAllstations, filterOrder='', metadata=['lowpass','LP'])

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
