#!/usr/bin/env python

##
## Grab the META data PKL generated from GetObsStations: (metapkl)
## Grab the Merged ADC,OBS,ERR time series data CSV computed by (mergedf) 
## Specify list up to four stations for comparison
##
import os,sys
import numpy as np
import pandas as pd
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates 
import seaborn as sns
#from utilities.utilities import utilities

##
## Build color/linetype mapping GLOBALS for known possible data types
##

colordict=dict()
colordict['Forecast']='b'
colordict['Nowcast']='b'
colordict['NOAA NOS']='g'
colordict['NOAA Tidal']='orange'
colordict['Difference']='r'
colordict['Contrails']='b'
colordict['NDBC']='g'
colordict['SWAN Forecast']='b'
colordict['SWAN Nowcast']='b'

dashdict=dict()
dashdict['Forecast']=(3,1)
dashdict['Nowcast']=(1,0)
dashdict['NOAA NOS']=(1,0)
dashdict['NOAA Tidal']=(1,0)
dashdict['Difference']=(1,0)
dashdict['Contrails']=(1,1)
dashdict['NDBC']=(2,2)
dashdict['SWAN Forecast']=(3,1)
dashdict['SWAN Nowcast']=(1,0)
dashdict['misc']=(3,2)

##
## Define mappings for y-axis plot titles
## Build more informative plots by better annotating the Legends
##

LEGENDS_MAP={
     'NOAA Tidal': 'NOAA Tidal: Water level',
     'NOAA NOS': 'NOAA NOS: Water level',
     'NDBC': 'NDBC: Wave height (WVHT)',
     'SWAN': 'SWAN: Wave height',
     'SWAN Forecast': 'Swan Forecast: Wave height',
     'SWAN Nowcast': 'Swan Nowcast: Wave height',
     'Contrails': 'Contrails: Water elevation (Stage)',
     'Forecast': 'Forecast: Water level',
     'Difference': 'Difference: Water level',
     'Nowcast': 'Nowcast: Water level'}

def union_all_source_stations(outputs_dict, station_id_list=None)->list:
    """
    The input dictionary contains a list of dataframes that are times x stations for the given source
    There are many opportunities for a specific station to be excluded from a given source such as a missing guage or a down server
    of missingness thresholding. Therefore, these dataframe may have different lists of stations.
    Here we take the union of all stations to simplify subsequent plotting

    If station_id_list is not None, then the returned data set is winnowed to only those 
    stations found in the station_id_list

    Parameters:
        outputs_dict:  a dictionary of times x stations objects with the currently known keys of:
        Forecast, Nowcast, NOAA NOS, NOAA Tidal, Difference, NDBC, Contrails 
        station_id_list: List (str) of stationids to retain

    Return:
        union dict of unique stations (str)
    """
    stations=list()
    for key,item in outputs_dict.items():
        #print(key)
        stations.extend(item.columns)
    station_list = list(set(stations))
    print('Union list of unique stations across all source objects is of length {}'.format(len(station_list)))
    if station_id_list is not None:
        station_list = { station for station in station_list if station in station_id_list }
        print('union_all_source_stations: Apply station_id_list to found station_list')
    return station_list

def map_datatype_to_plotter_type(in_key_name) -> str:
    """
    """ 
    station_types=['NOAA Tidal','NOAA NOS','Nowcast','Forecast', 'Difference']
    river_types=['Contrails']
    buoy_types=['NDBC','SWAN Forecast','SWAN Nowcast']

    if in_key_name in station_types:
        out_type='WATER_LEVEL'
    elif in_key_name in river_types:
        out_type='WATER_ELEVATION'
    elif in_key_name in buoy_types:
        out_type='WAVE_HEIGHT'
    else:
        out_type='None'
    return out_type

def union_station_types( outputs_dict, station_id_list=None) ->dict:
    """
    The input dictionary is used to try and classify the suite of per-station data
    into one of a few classes. These classes are used aftter-the-fact, by apsviz2 to 
    better inform the viewer if an indicated station is a station, buoy or river gauge
    The following procedure ASSUMES that no type overlapping occurs (eg, no ADCIRC WL with 
    an NDBC buoy output. However, not much checking is (or can) occur.

    Parameters:
        outputs_dict:  a dictionary of stations x type objects with the currently known keys of:
        Forecast, Nowcast, NOAA NOS, NOAA Tidal, Difference, NDBC, Contrails 
        station_id_list: List (str) of stationids to retain

    Return:
        union dict of stations and types (str). Defaults to the string: "None"
    """
##
## For each Key, fetch all the associated stations from the column list and appropriately assign the type
## Them check each LIST to see if more than one entry (besides nan) exist
##
    station_types=dict()
    for key,item in outputs_dict.items():
        datatype = map_datatype_to_plotter_type(key)
        dlist = item.columns.to_list() # bld a new dict with columns as the index
        for station in dlist:
            if station in station_types.keys():
                # Check and compare: Stop if inconsistency
                testtype=station_types[station]
                if testtype is 'None':
                    station_types[station]=datatype
                if testtype is not datatype:
                    utilities.log.error('Unexpected mixed station types found {} vs {} at station {}'.format(testtype,station_types[station],station))
                    sys.exit(1)
            else:
                station_types[station]=datatype
    return station_types

def union_station_names( outputs_meta_dict, station_id_list=None) ->dict:
    """
    Starting with a complete list of possible stations (derived from union_all_source_stations)
    compute the associated list of actual station name. These are used as Title for subsequent
    per-station plots. We only read OBSERVATIONS type data here. ADCIRC "stations" do not have a real name

    Parameters:
        outputs_meta_dict:  a dictionary of stations x metadata objects with the currently known keys of:
        LAT, LON, NAME, UNITS, TZ, OWNER, STATE, COUNTY

    Return:
        union list of unique stations (str) Optionally winnowed to the stations provided by station_id_list
        dict {stationid:stationname}
    """
    #possible_source_keys=['NOAA NOS','NOAA Tidal','Contrails','Forecast','Nowcast','SWAN FOrecast','SWAN Nowcast', 'NDBC'] # Preferred resolution order
    
    s_metadata=dict()
    for key in outputs_meta_dict.keys(): # possible_source_keys:
        #print(key)
        try:
            df_meta=outputs_meta_dict[key]
            for station in df_meta.index:
                df_meta.loc[station]
                lon,lat,name,state=df_meta.loc[station][['LON','LAT','NAME','STATE']]
                s_metadata[station]={'LON':lon,'LAT':lat,'NAME': name, 'STATE':state}
        except KeyError:
            pass
    if len(s_metadata)==0:
        print('No station names found at all {}. Abort'.format(outputs_meta_dict))
        sys.exit(1)
    print('Total number of station names found is {}'.format(len(s_metadata)))
    if station_id_list is not None:
        station_metadata = { station: s_metadata[station] for station in s_metadata.keys() if station in station_id_list }
    else:
        station_metadata = s_metadata
    return station_metadata

def get_bounds(df)->tuple:
    """
    The input dataset is expected to be for a single station and be a concatenation
    of all data sources that contributed. Here we try to find a good PRODUCT range for
    plotting. Most sources are periodic and range approx about zero. So we want a symmetric range
    Looks into every source to find their global min and global max

    We never expect NOAA/Contrails/ADCIRC water_levels (vs MSL) to be mixed with NDBC/SWAN
    wave_heights (vs 0). So try to adjust ymin for NDBC/SWAN to be zero

    Parameters
        df: time x station-sources (eg Nowcast, Forecast, etc)

    Returns:
        tuple of best ymin,ymax for this selected data collection
    """
    source_list = df.columns.tolist()
    ymin = abs(df.min().min())
    ymax = abs(df.max().max())

    if 'NDBC' in source_list or 'Contrails' in source_list or 'SWAN' in source_list:
        ymin=0.0
        ymax = ymax if ymax > ymin else -ymin
        ymax = (1.0+0.15)*ymax 
    else:
        # Look for the widest symmetric range if water_leels
        val = ymax if ymax > ymin else ymin
        ymin = -math.ceil(val)
        ymax = -ymin
    return(ymin, ymax)

def source_to_style(source_list):
    """
    Input an ordered source_list. The is generally, the order of COLUMNS of a 
    concat'd station-specifici dataframe with multiple sources. Here, we generate
    an ordered list of collor/line-style mappings for each source based on the GLOBAL:
    colordict, and dashdict

    Parameters
        source_list: list(str) of sources: Currently specified to be 
            Forecast, Nowcast, NOAA NOS, NOAA Tidal, Difference, Contrails 
    Return:
        colorlist: list(str), conforming order of color-to-source assignments
        dashlist: list(str), conforming order of linestyle-to-source assignments
    """
    colorlist = list()
    list_sources = list(colordict.keys())
    for source in source_list:
        try:
            color = colordict[source] 
        except KeyError:
            # No such source. Choose a generic color
            color = 'black'
        colorlist.append(color) 
    dashlist = list()
    list_sources = list(dashdict.keys())
    for source in source_list:
        try:
            dash = dashdict[source]
        except KeyError:
            # No such source. Choose a generic color
            print('No source {}'.format(source))
            dash = (2,2)
        dashlist.append(dash)
    return colorlist, dashlist

# Build a concatenated dataframe with headers of the indicated source

def build_source_concat_dataframe(outputs_dict, stationid)->pd.DataFrame:
    """
    For the input dictionary, grab each SOURCE and find the data for the desired station (if any)
    Merge results into a new dataframe with headers corresponding to the SOURCE

    NOTE: if one or more sources do not have the requested station those sources are summarily removed
    from the resulting df_station_specific_concat. I have observed cases where the NOAA obs does not have a station
    but the NOAA Tidal does (station=8770570). 

    Parameters
        stationid: (str) station id to extract data
        outputs_dict:  a dictionary of times x stations objects with the currewntly known keys of:
        Forecast, Nowcast, NOAA NOS, NOAA Tidal, Difference, Contrails

    Return
        station-specific dataframe (time x source)
    """
    list_dfs = list()
    for key,df in outputs_dict.items():
        try:
            list_dfs.append(df[stationid].to_frame(name=key))
        except KeyError:
            # No such station
            print('build_source_concat_dataframe: {}, no station_id {}: skip'.format(key,stationid))
    df_station_specific_concat = pd.concat(list_dfs,axis=1)
    return df_station_specific_concat

def update_headers_concat_dataframe(df_station_specific_concat)->pd.DataFrame:
    """
    For the input station-specific, source concatenated dataframe, replace the header SOURCE name 
    associated value from the global dict LEGENDS_MAP. These headers are used by the plotter routine
    to build legends for each plot. Since, different data asources may have different meanings (eg 
    water_level +/- max, vs water_height > 0), we needed that information reflected in the legend

    Headers are modified in-place

    Parameters:
        df_station_specific_concat: dataframe. time x source
    """
    header_list = df_station_specific_concat.columns.tolist()
    new_headers = [ LEGENDS_MAP[header] for header in header_list]
    df_station_specific_concat.columns = new_headers
    return df_station_specific_concat

# Cannot make any assumptons about the key order in the dict
def build_filename_map_to_csv(png_dict)->pd.DataFrame:
    """
    For APSVIZ convenience, instead of returning the station,lon,lat,name,filename,
    to the caller as a dictionary, return it as a dataframe.
    The intent being to save the data as a CSV instead of a JSON
    Build the converter code here, because this localizes the object key definitions

    Paramaters
        png_dict: A dict containing (at least) (LAT,LON,STATE,STATIONNAME,FIG)
    Return
        An equivelent dataframe with some data formatting checks
            Also, name case has changed slightly for backward compatibility
    """
    df=pd.DataFrame(png_dict['STATIONS']).T
    #df.columns=('Lat','Lon','State','StationName','Filename')
    df.rename(columns = {'LON':'Lon', 'LAT':'Lat', 'FILENAME':'Filename', 'TYPE':'Type', 'STATE':'State','STATIONNAME':'StationName'}, inplace = True)
    df.index.name='StationId'
    df['Lat'] = df['Lat'].astype(float)
    df['Lon'] = df['Lon'].astype(float)
    return df[['StationName','State','Lat','Lon','Filename','Type']]

# Build a station-specific plot
def create_station_specific_plot(fig, stationid, station_name, df_concat, time_range):
    """
    the input fig is updated in-place
    """
    tmin,tmax = time_range[0],time_range[1]
    ymin, ymax = get_bounds(df_concat)
    source_list = df_concat.columns.to_list()
    print(source_list)
    print(time_range)
    colorlist, dashlist = source_to_style(source_list)
    # Adjust headers for more informative plot legends
    df_concat = update_headers_concat_dataframe(df_concat)
    #
    sns.set_style('darkgrid')
    ax=sns.lineplot(data=df_concat, palette=colorlist, dashes=dashlist)
    ax.legend(loc = 4,fontsize = 6)
    ax.set_ylabel('meters', fontsize=7)
    ax.set_ylim([ymin, ymax])
    ax.set_xlim([tmin,tmax])
    ax.get_xaxis().set_visible(True)
    ax.set(xlabel=None)
    ax.xaxis.label.set_size(7)
    ax.tick_params(axis='x', labelsize=6)
    ax.tick_params(axis='y', labelsize=6)
    ax.grid(linestyle='-', linewidth='0.5', color='gray')
    plt.setp(ax.get_xticklabels(), rotation = 15)
    fig.suptitle(station_name)

def find_widest_timerange(outputs_dict):
    """
    This method identifies the global minimum and maximum times on the dataframe regardless of missingness
    Then, this range is passed to the plotter routine. This way all plots have the same width (x length)
    
    It is possible for a specific df to be completely empty and have NO time indexes. Especially, eg, 
        if grabbing old NDBC data. So we add a trap here to catch this scenario. Passing the all-nans data 
        forward is okay since the plotting routine will filter them out
    Parameter
        outputs_dict: dict of dfs, one for each source

    Return tmin,tmax in datetime64 format
    """
    time_mins = list()
    time_maxes = list()
    for key,df in outputs_dict.items():
        try:
            time_mins.append(min(df.index))
            time_maxes.append(max(df.index))
        except ValueError:
            # This would usually only happen to a completely empty dfs. Just skip and move on
            print('find_widest_timerange: Contains dfs with no data. Skip time range acquisition: Key ={}'.format(key))
            pass
    tmin = min(time_mins)
    tmax = max(time_maxes)
    return tmin,tmax

def generate_station_specific_PNGs(outputs_dict, outputs_meta_dict, outputdir='.', station_id_list=None )->pd.DataFrame:
    """
    The input dicts will be processed to build individual FIGS for each fund station id 
    Once a figure is created it cannot be put into a dict and then closed. So we MUST write out the file here.
    Let the calling program choose the output directory.

    Return
        figures_dict: A dictionary with the keys LAT,LON,STATE,STATIONNAME,FIG 
    """

    # Now grab the list of variables (eg ADC, OBS etc). Only need to choose a single station
    # How to handle DUPLICATE variable names?

    stations_union_source_list = union_all_source_stations(outputs_dict, station_id_list = station_id_list)
    station_names = union_station_names( outputs_meta_dict, station_id_list=station_id_list)
    station_types = union_station_types( outputs_dict, station_id_list = station_id_list)

    tmin,tmax = find_widest_timerange(outputs_dict)
    print('TIME RANGE {},{}'.format(tmin,tmax))

    data_dict=dict()
    for station in station_names.keys(): # Only want stations that have known ID names. 
        df_concat = build_source_concat_dataframe(outputs_dict, station)
        # Remove completely nan remaining stations
        df_concat.dropna(axis=1, how='all', inplace=True)
        print('After NAN reduction remaining sources is {}.Shape was {}'.format(station, df_concat.shape[1]))
        station_name = station_names[station]['NAME']
        lon = station_names[station]['LON']
        lat = station_names[station]['LAT']
        state = station_names[station]['STATE']
        st_type = station_types[station]
        # Build fig
        plt.close()
        fig = plt.figure(figsize=(6, 4))
        create_station_specific_plot(fig, station, station_name, df_concat, (tmin,tmax) )
        plt.subplots_adjust(bottom=0.25)
        filename=outputdir+'/'+station+'_WL.png'
        try:
            plt.savefig(filename)
            print('Wrote PNG to {}'.format(filename))
            data_dict[station]={'LAT':str(lat), 'LON':str(lon), 'STATE': state, 'STATIONNAME': station_name, 'TYPE': st_type, 'FILENAME':filename}
        except Exception as e:
            print('Failed writing PNG {}'.format(e))
        figures_dict = {'STATIONS':data_dict} # Conform to current apsviz2 procedure
        figure_dataframe = build_filename_map_to_csv(figures_dict)
    return figure_dataframe 
