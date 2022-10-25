#!/usr/bin/env python

import sys, os
import numpy as np
import pandas as pd
from utilities.utilities import utilities
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save_plot_model( adc_plot_grid=None, df_stations = None, df_surface=None, df_land_control=None, df_water_control=None, filename=None, plot_now=False, annotation=None):
    """
    Wrapper to take generated plt file and save to disk
    """
    if filename is None:
        utilities.log.error('save_plot_model: No filename provided to save the plot')
        return
    plt_data = plot_model(adc_plot_grid=adc_plot_grid, df_surface=df_surface, df_stations = df_stations, df_land_control=df_land_control, df_water_control=df_water_control, plot_now=plot_now, annotation=annotation)
    plt_data.savefig(filename, bbox_inches='tight')

def plot_model(adc_plot_grid=None, df_surface=None, df_stations=None, df_land_control=None, df_water_control=None, plot_now=True, annotation=None):
    """
    Basic plotter to display a 2D extrapolation field. The best images will include
    not only the surface, but also the stations, and control points as a reference 

    adc_plot_grid carries the (x,y) AXIS values spanning the data: df_surface
    So if, len(LAT) = 500, len(LON)=400, len(dim(df_surface)) is (400*500)

    Parameters:
        adc_plot_grid: dict with Keys of 'LON','LAT'
        df_surface: (DataFrame) (optional) with headers 'LON','LAT','VAL'
        df_land_control: (DataFrame) (optional) with headers 'LON','LAT','VAL'
        df_water_control: (DataFrame) (optional) with headers 'LON','LAT','VAL'
        plot_now: bool. True display the plot. Else not
    Results:
        A plot in the USA East Coast region
    """
    coastline=np.loadtxt('/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/coarse_us_coast.dat')
    cmap=plt.cm.jet
    print('annotate is {}'.format(annotation))
    #
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6,10), dpi=144) #, sharex=True)
    # Set up surface
    if df_surface is not None:
        print('plot_model: Found a extrapolated surface data set')
        x = adc_plot_grid['LON']
        y = adc_plot_grid['LAT']
        v = df_surface['VAL'].values
        v = v.reshape(-1,len(x)) # Reshapes to LAT column-major format
        mesh = ax.pcolormesh(x, y, v, shading='nearest', cmap=cmap, vmin=-0.2, vmax=0.2)
    if df_stations is not None:
        print('plot_model: Found a station data set')
        stations_X=df_stations['LON'].values
        stations_Y=df_stations['LAT'].values
        stations_V=df_stations['VAL'].values
        ax.scatter(stations_X, stations_Y, s=60, marker='o',
                   c=stations_V, cmap=cmap, edgecolor='black',
                   vmin=-0.2, vmax=0.2) # , cmap=cmap)
    # Merge control points
    if df_land_control is not None:
        print('plot_model: Found a land_control data set')
        land_X=df_land_control['LON'].values
        land_Y=df_land_control['LAT'].values
        land_V=df_land_control['VAL'].values
        ax.scatter(land_X, land_Y, s=30, marker='o',
                   c=land_V, cmap=cmap, edgecolor='black',
                   vmin=-0.2, vmax=0.2)
    if df_water_control is not None:
        print('plot_model: Found a water_control data set')
        utilities.log.info(f' water controls {df_water_control}')
        water_X=df_water_control['LON'].values
        water_Y=df_water_control['LAT'].values
        water_V=df_water_control['VAL'].values
        ax.scatter(water_X, water_Y, s=30, marker='x',
                   c='black', edgecolor='black') 
                   #, vmin=-0.2, vmax=0.2)
    ax.plot(coastline[:,0],coastline[:,1],color='black',linewidth=.25)
    ax.set_xlim([-90, -70])
    ax.set_ylim([20, 45])
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')
    plt.colorbar(mesh,orientation='vertical' )
    plt.title(annotation)
    if (plot_now):
        plt.show()
    return plt

def main(args):
    # Grab the test data
    f_land_control = '/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_land_controls.csv'
    f_water_control = '/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_water_controls.csv'
    f_surface = '/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_surface.csv'
    f_grid = '/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/adc_plot_grid.json'
     
    df_land_control=pd.read_csv(f_land_control,header=0,index_col=0)
    df_water_control=pd.read_csv(f_water_control,header=0, index_col=0)
    df_surface=pd.read_csv(f_surface,header=0, index_col=0)
    adc_coords = utilities.read_json_file(f_grid)
    plot_model(adc_plot_grid=adc_coords, df_surface=df_surface ,df_land_control=None,df_water_control=None)
    #plot_model(adc_plot_grid=adc_coords, df_surface=df_surface ,df_land_control=df_land_control,df_water_control=df_water_control)
    plt.show()

