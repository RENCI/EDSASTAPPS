#!/usr/bin/env python

##
## This modification of the compute_annual_errors.py simply reads OBS data that has been pregenerated/curated manually
## So the calls to NOAA have been removed. Everything else is basically the same
##
## Several data files are stored and in some cases, the same data are store in different formats (eg pkl and json). The pkl files are the 
## preferred format for working with subsequent python processors. But, the json formatted files facilitate users doing post procesing in 
## Matlab.   
##
## Select output path/filename are stored into a dict and written out to $RUNTIMEDIR/yearly_file_properties.json
##

import os,sys
import shutil
import numpy as np
import pandas as pd
import datetime as dt
import json
import re
import time as tm

import harvester.fetch_adcirc_data as fetch_adcirc_data
import harvester.generate_urls_from_times as genurls # generate_urls_from_times
import harvester.get_adcirc_stations as get_adcirc_stations
import harvester.get_observations_stations as get_obs_stations
import processing.compute_error_field as compute_error_field
import processing.interpolate_scaled_offset_field as interpolate_scaled_offset_field
import gridmap.grid_to_station_maps as grid_to_station_maps
import station_plotter as station_plotter

from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities
from argparse import ArgumentParser
# Fetch the station and Node data
# Must read our own pregenerated and detided obs data 

def remove_knockout_stations(knockout_dict, df_station) -> pd.DataFrame:
    """
    Occasionally, a station can report poor results for large ranges of times. This
    method allows the caller to nan those results
    Input should be a dataframe of time indexing and stations as columns.
    The time ranges specified in the args.knockout will be set to Nans inclusively.
    This method can be useful when the indicated statin has historically shown 
    poor unexplained performance over a large window of time

    Parameters:
       dataframe: time x stations.
       knockout_dict: <dict> 
    Returns:
       dataframe (time x stations) with knockout stations values set to nans
    """
    stations=list(knockout_dict.keys()) # How many. We anticipate only one but multiple stations can be handled
    cols=list(map(str,df_station.columns.to_list()))
    # 1 Do any stations exist? If not quietly leave
    if not bool(set(stations).intersection(cols)):
        return df_station
    # 2 Okay for each station loop over time range(s) and NaN away
    utilities.log.debug(f'Stations is {stations}')
    utilities.log.debug(f'dict {knockout_dict}')
    print(df_station)
    for station in stations:
        for key, value in knockout_dict[station].items():
            df_station.loc[value[0]:value[1]][station]=np.nan
    utilities.log.info('Knockout dict has been applied')
    return df_station


def get_adcirc_stations_fort63_style(fname=None, station_filename='noaa_station_file'):
    if fname is None:
        utilities.log.error('No Fort63_style ADCIRC station file assigned: Abort')
        sys.exit(1)
    df = pd.read_csv(fname, index_col=0, header=0, skiprows=[1], sep=',',dtype=str)
    print(df)
    try:
        df["Node"]=df["Node"].astype(int)
    except IndexError as e:
        utilities.log.error(f'Unsuccessful fort63_station read. Perhaps no Node column ? {e}')
    return df

    if args.fort63_style:
        adcirc_stations=get_adcirc_stations_fort63_style(fname=station_filename)
    else:
        adcirc_stations=get_adcirc_stations_fort61_style(fname=station_filename)
    adcirc_metadata='Raw_data'
    print(adcirc_stations)
    print(adcirc_metadata)
    return adcirc_stations

# noinspection PyPep8Naming,DuplicatedCode
def main(args):
    t0 = tm.time()

    #assumes the existance of a proper main.yml to get IO information
    if args.main_yamlname is None:
        main_yamlname=os.path.join(os.path.dirname(__file__), './config', 'main.yml')
    else:
        main_yamlname=args.main_yamlname
    config = utilities.init_logging(subdir=None, config_file=main_yamlname, log_file_metadata='annual')

    if args.obs_filename is None:
        utilities.log.error('Must specify an input OBS data file')
        sys.exit(1)
    obs_filename=args.obs_filename
    utilities.log.info(f'OBS input filename is {obs_filename}')

    if args.meta_filename is None:
        utilities.log.error('Must specify an input meta data file')
        sys.exit(1)
    meta_filename=args.meta_filename
    utilities.log.info(f'META input filename is {meta_filename}')

    #
    grid = args.gridname 
#
# Specify the details needed for processing the URL as a "raw" (aka local) url
# We need to be prepared to pass a NONE grid to AST to get it to read the URL
# without tryinmg to parse the grid/instance/site names from it
#
    print('Raw url data structure expected')
    utilities.log.info('Running a RAW url type')
    ensemble = 'NONE'
    gridname = 'NONE' # This is intentionally set to NONE to trigger Raw url lookup
    sitename = 'NONE'

    variable_name = args.variable_name
    utilities.log.info('Chosen variable_name is {}'.format(variable_name))

    if args.fort63_style:
        utilities.log.info('Fort_63 style station inputs specified')

    if args.knockout is not None:
        knockout_dict = io_utilities.read_json_file(args.knockout)
        utilities.log.info('Read a knockout file')
    else:
        knockout_dict=None

    # fetch the grid-specific station data. Land/water control points not (yet) needed for apsviz2 work
    station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.gridname, mapfile=args.map_file, datatype='NOAA_STATIONS')
    file_land_controls = grid_to_station_maps.find_land_control_points_from_map(gridname=args.gridname, mapfile=args.map_file)
    file_water_controls = grid_to_station_maps.find_water_control_points_from_map(gridname=args.gridname, mapfile=args.map_file)

    print(f' station file {station_file}')

    # Ensure the input URL is of the proper style by converting it as necessary (63 vs 61)
    input_url = args.url
    print(f' input url {input_url}')
    if args.fort63_style:
        if args.custom_fort63_name is None:
            input_url=get_adcirc_stations.convert_urls_to_63style([args.url])[0] # Assumed to only have one of them
        else:
            input_url=get_adcirc_stations.convert_urls_to_63style_customfilename([args.url], filename=args.custom_fort63_name)[0]
            utilities.log.info('Using a customized name for the fort.63.nc filename')
    else:
        input_url=get_adcirc_stations.convert_urls_to_61style([args.url])[0]
    print(input_url)

##
# Get the station file name
# Skip the original method. We assume the data are fort63_compilant and state so
# Need to do this differenty
# noaa_station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.gridname, mapfile=args.map_file, datatype='NOAA_STATIONS')

    print(f'meta filename {meta_filename}')
    print(meta_filename)
    adcirc_stations = get_adcirc_stations_fort63_style(fname=meta_filename)
    print(f' ADCIRC stations {adcirc_stations}')

#
# Holders for plotting data
#
    outputs_dict=dict()
    outputs_metadict=dict()
#
# Fetch the actual ADCIRC data. No need to use knockout here since compute_error will intersect times anyway
#
    df_adcirc_data, df_adcirc_meta = fetch_adcirc_data.process_adcirc_stations( [input_url], adcirc_stations, gridname, ensemble, sitename, 'water_level', resample_mins=60, fort63_style=args.fort63_style, variable_name=variable_name)
    df_adcirc_data.replace('-99999',np.nan,inplace=True)
    df_adcirc_meta.replace('-99999',np.nan,inplace=True)

    outputs_dict['Nowcast']=df_adcirc_data
    outputs_metadict['Nowcast']=df_adcirc_meta
    utilities.log.info('Finished with ADCIRC')

# What YEAR does this data apply to? Find it and then build a time range for the observations

    timein = min(df_adcirc_data.index)
    timeout = max(df_adcirc_data.index)

    starttime = timein.strftime('%Y-%m-%d %H:%M:%S')
    endtime = timeout.strftime('%Y-%m-%d %H:%M:%S')
    time_range=(starttime,endtime)
#
# Fetch the observations data
# Detailed data is collected at maximum frequency. User resample is applied to subsequent smoothing
#
    obs_filename=args.obs_filename
    print(f' Try to open {obs_filename}')
    #inobsdir='/projects/sequence_analysis/vol1/prediction_work/NOAA.Reanalysis/43_OBSERVATIONS_BRIAN_DETRENDING'
    #inobsfile='noaa_obs_1979_2022_detrended'

    df_all_obs_data_full=pd.read_csv(obs_filename, header=0, index_col=0)
    df_all_obs_data_full.index = pd.to_datetime(df_all_obs_data_full.index)
    
#NEED TO KNOCKOUT
    if knockout_dict is not None:
        print('Applying a knockout to the OBS dataframe')
        df_all_obs_data_full = remove_knockout_stations(knockout_dict, df_all_obs_data_full.copy())


    ##############################################################
    # This is not right
    # What year are in looking at?

    year = timein.year
    print(f'Year determined from ADCIRC to be {year}')
    df_all_obs_data = df_all_obs_data_full.loc[str(year)]

## Must change this kludged in bcs of new,local OBS data source
    df_all_meta = pd.read_csv(meta_filename, index_col=0, header=0, skiprows=[1],dtype=str)
    df_all_meta['lon']=df_all_meta['lon'].astype('float')
    df_all_meta['lat']=df_all_meta['lat'].astype('float')
    df_all_meta=df_all_meta[['stationid','lon','lat','stationname']].rename(columns={'stationid': 'STATION', 'lon': 'LON', 'lat':'LAT', 'stationname':'NAME' })
    df_all_meta.set_index('STATION',inplace=True)
    meta_obs = df_all_meta.copy()

    print(df_all_obs_data)
    print(df_all_meta)

    meta_obs_list=df_all_obs_data.columns.tolist()

    utilities.log.info(f' STATIC OBS list {df_all_obs_data.columns.tolist()}')
    utilities.log.info(f' STATIC META list {meta_obs.index.to_list()}')

    meta_thresholded = meta_obs.loc[meta_obs_list]

#
# Apply a moving average (smooth) the data performed the required resampling to the desired rate followed by interpolating
# Note for hourly data we are imposing a window of 11 hours wide... centered
#

    df_observations_data = df_all_obs_data
    df_observations_meta = meta_thresholded

    outputs_dict['NOAA NOS']=df_observations_data # df_observations_data
    outputs_metadict['NOAA NOS']=df_observations_meta
    utilities.log.info('Finished with NOAA Observations')

#
# Fetch the intersection of stations
#
    data_station_set = set(df_adcirc_data.columns.tolist()).intersection(df_observations_data.columns.tolist())
    meta_station_set = set(df_adcirc_meta.index.tolist()).intersection(df_observations_meta.index.tolist())
    final_station_list = list(set(data_station_set).intersection(meta_station_set))

#
# Compute errors for all stations
#
    comp_err = compute_error_field.compute_error_field(df_observations_data, df_adcirc_data, df_observations_meta, n_hours_per_period=24) 
    comp_err._intersection_stations()
    comp_err._intersection_times()
    comp_err._compute_and_average_errors()
    # Remove outliers ?
    #utilities.log.info('Attempt to remove outliers')
    #comp_err._remove_station_outliers()
    utilities.log.info('No Z-score outlier removal')

    outputs_dict['Difference']=comp_err.diff
    utilities.log.info('Finished with Compute Errors')
    print(f' DIFF {comp_err.diff}')

# Generate a set of PNGs from the plotter code

    rootdir=io_utilities.construct_base_rootdir(config['DEFAULT']['RDIR'], base_dir_extra=None)
    utilities.log.info('Output directory specified to be {}'.format(rootdir))

#
# Generate PER_STATION plots for APSVIZ2 insets
# This is failing
#
    #print(outputs_dict)
    #print(outputs_metadict)
    #df_station_file_png_locations = station_plotter.generate_station_specific_PNGs(outputs_dict,
    #    outputs_metadict, outputdir=rootdir, station_id_list=None)
    #print(f'Locations {df_station_file_png_locations}')
    #station_props = io_utilities.write_csv(df_station_file_png_locations, rootdir=rootdir,subdir='',fileroot='stationProps',iometadata='')
    #utilities.log.info('Wrote out station_properties file to {}'.format(station_props))

#
# Save data files for daily processing and interpolation
# Build a run_file_properties.dict that point to all fully qualified files names generated by this method
# By default store it into $RUNTIMEDIR
#

    output_files=dict()

# Save ADCIRC fullgrid coordinates

    if args.custom_fort63_name is None:
        url_63=get_adcirc_stations.convert_urls_to_63style([input_url])[0] # Assumed to only have one of them
    else:
        url_63=get_adcirc_stations.convert_urls_to_63style_customfilename([input_url], filename=args.custom_fort63_name)[0]
        utilities.log.info('Using a customized name for the fort.63.nc filename')
    #print(url_63)
    iosubdir='' # 'adcpkl'
    iometadata=''
    adc_coords = get_adcirc_stations.extract_adcirc_grid_coords(url_63 )
    ADCfilecoordsJson = io_utilities.write_dict_to_json(adc_coords, rootdir=rootdir,subdir=iosubdir,fileroot='adc_coord',iometadata=iometadata)
    io_utilities.write_json_file(adc_coords, ADCfilecoordsJson)
    utilities.log.info('Wrote grid coords to {}'.format(ADCfilecoordsJson))
    output_files['adc_coord']=ADCfilecoordsJson

# Data from the compute errors routines
    iosubdir='' # 'errorfield'
    iometadata=''
    # Write selected Pickle data 
    tideTimeErrors = io_utilities.write_pickle(comp_err.diff,rootdir=rootdir,subdir=iosubdir,fileroot='tideTimeErrors',iometadata=iometadata)
    utilities.log.info('Wrote out Pickle {}'.format(tideTimeErrors))
    output_files['tideTimeErrors']=tideTimeErrors

    # Write selected Pickle data 

    # Write selected in JSON format. Basically the Merged_dict data multi-index set
    data_dict = compute_error_field.combine_data_to_dict(comp_err.adc,comp_err.obs,comp_err.diff, product='WL')
    dict_json = io_utilities.write_dict_to_json(data_dict, rootdir=rootdir,subdir=iosubdir,fileroot='adc_obs_error_merged',iometadata=iometadata)
    utilities.log.info('Wrote out JSON {}'.format(dict_json))
    output_files['adc_obs_error_merged']=dict_json

    # Write the CSV files which carry averages
    station_summary_aves=io_utilities.write_csv(comp_err.df_final, rootdir=rootdir,subdir=iosubdir,fileroot='stationSummaryAves',iometadata=iometadata)
    station_period_aves=io_utilities.write_csv(comp_err.df_cycle_aves, rootdir=rootdir,subdir=iosubdir,fileroot='stationPeriodAves',iometadata=iometadata)
    utilities.log.info('Wrote out CSV cyclic averages {}'.format(station_period_aves))
    utilities.log.info('Wrote out CSV summary averages {}'.format(station_summary_aves))
    output_files['stationSummaryAves']=station_summary_aves
    output_files['stationPeriodAves']=station_period_aves

# Raw water level data OBSERVATIONS
    iosubdir='' # 'obspkl'
    iometadata='hourly_height'
    metapkl = io_utilities.write_pickle(df_observations_meta,rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_metadata',iometadata=iometadata)
    detailedpkl = io_utilities.write_pickle(df_observations_data, rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_detailed',iometadata=iometadata)
    output_files['obs_wl_detailed_pkl']=detailedpkl
    output_files['obs_wl_metadata_pkl']=metapkl

# Also write out in JSON format. Eg for use by our Matlab processors

    df_observations_data.index = df_observations_data.index.strftime('%Y-%m-%d %H:%M:%S')
    detailedjson = io_utilities.write_json(df_observations_data,rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_detailed',iometadata=iometadata)
    metajson = io_utilities.write_json(df_observations_meta,rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_metadata',iometadata=iometadata)
    output_files['obs_wl_detailed_json']=detailedjson
    output_files['obs_wl_metadata_json']=metajson

# Raw ADCIRC water level
    metapkl = io_utilities.write_pickle(df_adcirc_meta,rootdir=rootdir,subdir='',fileroot='adc_wl_metadata',iometadata='')
    detailedpkl = io_utilities.write_pickle(df_adcirc_data, rootdir=rootdir,subdir='',fileroot='adc_wl_detailed',iometadata='')
    output_files['adc_wl_detailed']=detailedpkl
    output_files['adc_wl_metadata']=metapkl

# Write out the dict file for subsequent processing
    output_files_json = io_utilities.write_dict_to_json(output_files, rootdir=rootdir,subdir='',fileroot='yearly_file_properties',iometadata='')
    utilities.log.info('Wrote out dict: {}'.format(output_files_json))
    utilities.log.info('Finished: Total runtime was {}s'.format(tm.time()-t0))

if __name__ == '__main__':
    from argparse import ArgumentParser
    import sys

    parser = ArgumentParser()

    parser.add_argument('--rootdir', action='store', dest='rootdir', default=None,
                        help='Available high leverl directory')
    parser.add_argument('--ignore_pkl', help="Ignore existing pickle files.", action='store_true')
    parser.add_argument('--doffset', default=-4, help='Day lag or datetime string for analysis: def to YML -4', type=int)
    parser.add_argument('--iometadata', action='store', dest='iometadata',default='', help='Used to further annotate output files', type=str)
    parser.add_argument('--iosubdir', action='store', dest='iosubdir',default='', help='Used to locate output files into subdir', type=str)
    parser.add_argument('--urljson', action='store', dest='urljson', default=None,
                        help='String: Filename with a json of urls to loop over.')
    parser.add_argument('--url', action='store', dest='url', default=None,
                        help='String: url.')
    parser.add_argument('--grid', default='hsofs',dest='grid', help='Choose name of available grid',type=str)
    parser.add_argument('--knockout', default=None, dest='knockout', help='knockout jsonfilename', type=str)
    parser.add_argument('--fort63_style', action='store_true',
                        help='Boolean: Will inform Harvester to use fort.63.methods to get station nodesids')
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--data_source', action='store', dest='data_source', default='TDS', type=str,
                        help='choose supported data source: default = TDS')
    parser.add_argument('--map_file', action='store',dest='map_file', default=None,
                        help='str: Select appropriate map_file yml for grid lookup')
    parser.add_argument('--main_yamlname', action='store',dest='main_yamlname', default=None,
                        help='str: Select appropriate main_yamlname')
    parser.add_argument('--gridname', action='store',dest='gridname', default='hsofs',
                        help='str: Select appropriate gridname Default is hsofs')
    parser.add_argument('--variable_name', action='store',dest='variable_name', default='zeta',
                        help='str: Select non default ADCIRC NC dataset variable name')
    parser.add_argument('--custom_fort63_name', action='store',dest='custom_fort63_name', default=None,
                        help='str: Select an alternative name for the fort.63.nc analog filename Usually when using rechunked data')
    parser.add_argument('--obs_filename', action='store',dest='obs_filename', default=None,
                        help='str: FQFN for the pregenerated OBS observational data')
    parser.add_argument('--meta_filename', action='store',dest='meta_filename', default=None,
                        help='str: FQFN for the pregenerated META observational data')

    args = parser.parse_args()
    sys.exit(main(args))
