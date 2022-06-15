#!/usr/bin/env python

import os,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities

def make_plot(start, end, station, src, stationName, dfs, dfs_7d, dfs_weekly_mean, dfs_monthly_mean, odir): 
    sns.set(rc={'figure.figsize':(11, 4)}) # Setr gray background and white gird
    col = sns.color_palette("bright")[0:3] # CHanged to make figs for reananalysis more consistent

    # Plot daily, weekly resampled, and 7-day rolling mean time series together
    fig, ax = plt.subplots()
    ax.plot(dfs.loc[start:end, src],
    marker='.', markersize=1, linewidth=0.1,label='Hourly', color=col[0]) # color='gray',)
    ax.plot(dfs_7d.loc[start:end, src],
    color='red', alpha=0.3, linewidth=.5, linestyle='-', label='7-d Rolling Mean') # red
    ax.plot(dfs_weekly_mean.loc[start:end, src],
    marker='o', color=col[2] ,markersize=6, linestyle='-', label='Weekly Mean') # green
    ax.plot(dfs_monthly_mean.loc[start:end, src],
    color='black',linewidth=0.5, linestyle='-', label='Monthly Mean') # black
    ax.set_ylabel('WL [m]')
    ax.set_title(stationName, fontdict={'fontsize': 12, 'fontweight': 'medium'})
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=None))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'));
    ax.legend(fontsize=10);
    plt.xticks(rotation=0, fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    fname='/'.join([odir,station+'.png'])
    #fname=station+'.png'
    utilities.log.info('Saving station png {}'.format(fname))
    plt.savefig(fname)
    plt.close()
    #plt.show()

# make_lowpass_plot no longer works
def make_lowpass_plot(start, end, lowpassAllstations, filterOrder='', metadata=['lowpass','LP']):
    """
    On entry the dict (lowpassAllStations) carries all the OBS,DIFF,LP data sets
    Plot the OBS and Diff data in a prominant way. Then layer on the cutoffs
    """
    col = sns.color_palette("bright")[0:3] 
    colMetadata = metadata[1]
    nameMetadata = metadata[0]
    plotterDownmove=0.6
    station=lowpassAllstations['station']
    stationName=lowpassAllstations['stationName']
    print('station {} Name {}'.format(station, stationName))
    #
    plt.close()
    sns.set(rc={'figure.figsize':(11, 4)}) # Set gray background and white gird
    fig, ax = plt.subplots()
    # OBS and DIFF
    ax.plot(lowpassAllstations[station]['OBS'][start:end],
    marker='.', markersize=1, linewidth=0.1,color=col[0],label='Obs Hourly') # gray
    ax.plot(lowpassAllstations[station]['ERR'][start:end],
    color=col[2], marker='o',markersize=2, linewidth=.5, linestyle='-', label='ADCIRC-OBS') # black
    # Now start with the lowpass filters
    df_lp = lowpassAllstations[station][colMetadata]
    cutoffs=df_lp.columns.to_list()
    shiftDown=0.0
    for cutoff in cutoffs:
        shiftDown+=plotterDownmove
        ax.plot(df_lp[cutoff][start:end]-shiftDown,
        marker='x', markersize=2, linewidth=0.1,label='_'.join([nameMetadata,cutoff]))
    #
    ax.set_ylabel('WL [m]')
    if filterOrder=='':
        ax.set_title(stationName+'. FFT', fontdict={'fontsize': 12, 'fontweight': 'medium'})
    else:
        ax.set_title(stationName+'. Polynomial Order='+str(filterOrder), fontdict={'fontsize': 12, 'fontweight': 'medium'})
    ax.xaxis.set_major_locator(mdates.MonthLocator(bymonth=None))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'));
    ax.legend(fontsize=10);
    plt.xticks(rotation=0, fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    fileMetaname = '_'.join([str(station),nameMetadata])
    fileName = fileMetaname+'.png' if filterOrder=='' else fileMetaname+'_'+str(filterOrder)+'.png'
    plt.savefig(fileName)
    #plt.show()
