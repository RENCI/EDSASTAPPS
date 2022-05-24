#!/usr/bin/env python

##
## This implementation of the apsviz2 code (formally) is implemented to use the AST infrastructure
## This is designed to only support the template_url harvester approach to building URLs
##
## While the code is resilient to missingness in nowcasts, and observations, the starting forecast
## must be available otherwise the code will hard-exit.
##

import os,sys
import shutil
import numpy as np
import pandas as pd
import datetime as dt

import harvester.fetch_adcirc_data as fetch_adcirc_data
import harvester.generate_urls_from_times as genurls # generate_urls_from_times
import harvester.get_adcirc_stations as get_adcirc_stations
import harvester.get_observations_stations as get_obs_stations
import processing.compute_error_field as compute_error_field
import processing.interpolate_scaled_offset_field as interpolate_scaled_offset_field
import station_plotter as station_plotter
import gridmap.grid_to_station_maps as grid_to_station_maps

from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities
from argparse import ArgumentParser

##
## Run the canonical APSVIZ2 procedure
##

dformat='%Y-%m-%d %H:%M:%S'

def main(args):
    """
    assumes the existance of a proper main.yml to get IO information
    """
    if args.main_yamlname is None:
        main_yamlname=os.path.join(os.path.dirname(__file__), './config', 'main.yml')
    else:
        main_yamlname=args.main_yamlname
    config = utilities.init_logging(subdir=None, config_file=main_yamlname)

    # Set up IO env
    utilities.log.info("Product Level Working in {}.".format(os.getcwd()))
    utilities.log.info('Selected main yaml file is {}'.format(main_yamlname))

    outputs_dict=dict()
    outputs_metadict=dict()

    # Basic checks
    if args.config_name is None:
        config_name =os.path.join(os.path.dirname(__file__), './secrets', 'url_framework.yml')

    # Extra bits incase an hurricane URL is passed in 
    hurricane_source=args.hurricane_source
    hurricane_year=args.hurricane_year

    # Generally we expect this to be False. sinbce many of the stations will not have ADCIRC Nodeid data
    fort63_style=args.fort63_style # Need to see if we actually want to use fort63 approaches

    # Specify stations (NOAA) to use for ADCIRC sources

    # fetch the grid-specific station data. Land/water control points not (yet) needed for apsviz2 work
    station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.gridname, mapfile=args.map_file, datatype='NOAA_STATIONS')
    #file_land_controls = grid_to_station_maps.find_land_control_points_from_map(gridname=args.gridname, mapfile=args.map_file)   # ,mapfile=map_file)
    #file_water_controls = grid_to_station_maps.find_water_control_points_from_map(gridname=args.gridname, mapfile=args.map_file)

    # Ensure the input URL is of the proper style (63 vs 61)
    input_url = args.url
    if fort63_style:
        input_url=get_adcirc_stations.convert_urls_to_63style([args.url])[0] # Assumed to only have one of them
    else:
        input_url=get_adcirc_stations.convert_urls_to_61style([args.url])[0]
    print(input_url)

##
## ADCIRC step 1. Process the INPUT (Forecast) URL. 60min default sampling
##
    adc = get_adcirc_stations.get_adcirc_stations(source='ASGS', product=args.data_product,
                station_list_file=station_file,
                knockout_file=None, fort63_style=fort63_style )
    data_adc,meta_adc=adc.fetch_station_product( [input_url] , return_sample_min=args.return_sample_min, fort63_style=fort63_style  )
    # Revert Harvester filling of nans to -99999 back to nans
    data_adc.replace('-99999',np.nan,inplace=True)
    meta_adc.replace('-99999',np.nan,inplace=True)
    outputs_dict['Forecast']=data_adc
    outputs_metadict['Forecast']=meta_adc
    utilities.log.info('Finished with ADCIRC Forecasts')

    # Grab the stop and start times from the data set. Will be needed for tidal predictions data
    time_index=data_adc.index.tolist()
    starttime = min(time_index).strftime('%Y-%m-%d %H:%M:%S')
    endtime = max(time_index).strftime('%Y-%m-%d %H:%M:%S')

##
## ADCIRC step 2. Process the associated NOWCAST URLs.
##

# starttime/endtime can be > resulting data set (urls list) if not enough ADCIRC files were actually found.
# We still want the full range for observations and tidal predictions

    # Need to build a set of NOWCASTs from the input url.
    nowadc = genurls.generate_urls_from_times(url=input_url,ndays=args.ndays)
    now_urls = nowadc.build_url_list_from_template_url_and_offset(ensemble='nowcast')
    print(now_urls)
    nowadc = get_adcirc_stations.get_adcirc_stations(source='ASGS', product=args.data_product,
                station_list_file=station_file,
                knockout_file=None, fort63_style=fort63_style )
    data_now_adc,meta_now_adc=nowadc.fetch_station_product(now_urls, return_sample_min=args.return_sample_min, fort63_style=fort63_style  )
    data_now_adc.replace('-99999',np.nan,inplace=True)
    meta_now_adc.replace('-99999',np.nan,inplace=True)
    outputs_dict['Nowcast']=data_now_adc
    outputs_metadict['Nowcast']=meta_now_adc
    utilities.log.info('Finished with ADCIRC Nowcasts')
    print(now_urls)

    # Figure out what the asociated range of times is. 
    obs_endtime = starttime # '%Y-%m-%d %H:%M:%S': Grab the beginning of the forecast
    dt_starttime = dt.datetime.strptime(obs_endtime,dformat)+dt.timedelta(days=args.ndays) # How many days BACK
    obs_starttime=dt.datetime.strftime(dt_starttime, dformat)
    print('Nowcast time range is from {} through {}'.format(obs_starttime, obs_endtime))

##
## Assemble the list of valid observations from which to compute error series'
##
    valid_obs=list()
    valid_obs_meta=list()

## OBSERVATIONS. Get the NOAA station data
##
# Detailed data is collected at maximum frequency. User resample is applied to subsequent smoothing

    noaa_station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.gridname, mapfile=args.map_file, datatype='NOAA_STATIONS')
    if noaa_station_file is not None:
        obs = get_obs_stations.get_obs_stations(source='NOAA', product='water_level',
                station_list_file=noaa_station_file)
        data_obs,meta_obs=obs.fetch_station_product((obs_starttime,obs_endtime), return_sample_min=0, interval='None' )
        data_obs.replace('-99999',np.nan,inplace=True)
        meta_obs.replace('-99999',np.nan,inplace=True)
        # Remove stations with too many nans ( Note Harvester would have previously removed stations that are ALL NANS)
        data_thresholded = obs.remove_missingness_stations(data_obs, max_nan_percentage_cutoff=10)  # (Maximum allowable nans %)
        meta_thresholded = meta_obs.loc[data_thresholded.columns.tolist()]
        meta_obs_list = set(data_thresholded.columns.tolist()).intersection(meta_obs.index.to_list())
        meta_obs_thresholded = meta_obs.loc[meta_obs_list]
        # Apply a moving average (smooth) the data performed the required resampling to the desired rate followed by interpolating
        data_obs_smoothed = obs.fetch_smoothed_station_product(data_thresholded, return_sample_min=args.return_sample_min, window=11)
        outputs_dict['NOAA NOS']=data_obs_smoothed
        outputs_metadict['NOAA NOS']=meta_obs_thresholded
        valid_obs.append(data_obs_smoothed)
        valid_obs_meta.append(meta_obs)
        utilities.log.info('Finished with NOAA Observations')

##
## OBSERVATIONS #2. Get the Contrails station data. Assumes the contrails stations are in the same station_list
## Only grabs RIVER data

    utilities.log.info('Preparing for a CONTRAILS fetch')
    contrails_config = args.contrails_auth 
    utilities.log.info('Got Contrails access information {}'.format(contrails_config))
    rivers_data_product='river_water_level'
    contrails_station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.gridname, mapfile=args.map_file, datatype='CONTRAILS_RIVERS')
    if contrails_station_file is not None:
        contrails = get_obs_stations.get_obs_stations(source='CONTRAILS', product=rivers_data_product,
            contrails_yamlname=contrails_config, station_list_file=contrails_station_file)
        # Get data at highest resolution. Return at 15min intervals
        cont_data,cont_meta=contrails.fetch_station_product((obs_starttime,obs_endtime), return_sample_min=args.return_sample_min, interval='None' )
        cont_data.replace('-99999',np.nan,inplace=True)
        cont_meta.replace('-99999',np.nan,inplace=True)
        cont_meta_list = set(cont_data.columns.tolist()).intersection(cont_meta.index.to_list())
        cont_meta = cont_meta.loc[cont_meta_list]
        outputs_dict['Contrails']=cont_data
        outputs_metadict['Contrails']=cont_meta
        valid_obs.append(cont_data)
        valid_obs_meta.append(cont_meta)
        utilities.log.info('Finished with Contrails Observations')
##
## OBSERVATIONS. #3 Get known Buoy wave_height data
##
# We expect a lot of missingness for these data so set a high number=90% Do not use 100% because if your entire dataset is empty
# (eg a wrong time range) then all empty buoys will be kept

    ndbc_station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.gridname, mapfile=args.map_file, datatype='NDBC_BUOYS')
    if ndbc_station_file is not None:
        ndbc = get_obs_stations.get_obs_stations(source='NDBC', product='wave_height',
                station_list_file=ndbc_station_file)
        data_ndbc,meta_ndbc=ndbc.fetch_station_product((obs_starttime,obs_endtime), return_sample_min=0, interval='None' )
        data_ndbc.replace('-99999',np.nan,inplace=True)
        meta_ndbc.replace('-99999',np.nan,inplace=True)
        # Remove stations with too many nans ( Note Harvester would have previously removed stations that are ALL NANS)
        try:
            data_ndbc_thresholded = ndbc.remove_missingness_stations(data_ndbc, max_nan_percentage_cutoff=90)  # (Maximum allowable nans %)
            meta_thresholded = meta_ndbc.loc[data_ndbc_thresholded.columns.tolist()]
            meta_ndbc_list = set(data_ndbc_thresholded.columns.tolist()).intersection(meta_ndbc.index.to_list())
            meta_ndbc_thresholded = meta_ndbc.loc[meta_ndbc_list]
            outputs_dict['NDBC']=data_ndbc_thresholded
            outputs_metadict['NDBC']=meta_ndbc_thresholded
            utilities.log.info('Finished with NDBC Observations')
        except Exception as e:
            utilities.log.error('NDBC: Failed to find any data: do not include to plot list {}'.format(e))
        valid_obs.append(data_ndbc)
        valid_obs_meta.append(meta_ndbc)
##
## Combine observations data for the error computations - No NDBC data here
##
    data_obs_smoothed = pd.concat(valid_obs,axis=1)
    meta_obs = pd.concat(valid_obs_meta,axis=0)
    print(meta_obs)

##
## TIDAL PREDICTIONS. Get the NOAA station data
##
    # Want times from beginning of the nowcasts to the end of the ADCIRC forecasts

    pred = get_obs_stations.get_obs_stations(source='NOAA', product='predictions',
                station_list_file=station_file)

    data_pred,meta_pred=pred.fetch_station_product((obs_starttime,endtime), return_sample_min=args.return_sample_min, interval='None' )
    data_pred.replace('-99999',np.nan,inplace=True)
    meta_pred.replace('-99999',np.nan,inplace=True)
    meta_pred_list = set(data_pred.columns.tolist()).intersection(meta_pred.index.to_list())
    meta_pred = meta_pred.loc[meta_pred_list]

    outputs_dict['NOAA Tidal']=data_pred
    outputs_metadict['NOAA Tidal']=meta_pred
    utilities.log.info('Finished with Tidal Predictions')

##
## Perform the error computations
##

    comp_err = compute_error_field.compute_error_field(data_obs_smoothed, data_now_adc, meta_obs) # All default params
    comp_err._intersection_stations()
    comp_err._intersection_times()
    comp_err._compute_and_average_errors()

    outputs_dict['Difference']=comp_err.diff
    utilities.log.info('Finished with Compute Errors')

##
## Select from RUNTIMEDIR, where to write these files
##

    print('POOP')
    print(config['DEFAULT']['RDIR'])
    print(type(config['DEFAULT']['RDIR']))
    print(args.finalDIR)

    if args.finalDIR is None:
        rootdir=io_utilities.construct_base_rootdir(config['DEFAULT']['RDIR'], base_dir_extra=None)
    else:
        print('Override with finalDIR setting {}'.format(args.finalDIR))
        rootdir=io_utilities.construct_base_rootdir(args.finalDIR, base_dir_extra=None)
    utilities.log.info('Output directory specified to be {}'.format(rootdir))

##
## Generate PER_STATION plots for APSVIZ2 insets
##
    df_station_file_png_locations = station_plotter.generate_station_specific_PNGs(outputs_dict, 
        outputs_metadict, outputdir=rootdir, station_id_list=None)
    print(df_station_file_png_locations)
    station_props = io_utilities.write_csv(df_station_file_png_locations, rootdir=rootdir,subdir='',fileroot='stationProps',iometadata='')
    utilities.log.info('Wrote out station_properties file to {}'.format(station_props))

##
## Copy over logfile into the rootdir so as to reside with the output data sets
##
    shutil.copy(utilities.LogFile,'/'.join([rootdir,'logs'])) # Copy and rename to logs for apsviz2 pipeline to find
    utilities.log.info('Copy log file')

##
## construct and save the png CSV lookup file for APSVIZ2
##
    # Save the data as test data
    #for key,item in outputs_dict.items():
    #    f=key.replace(' ','_')
    #    item.to_pickle(f+'.pkl')
    #    print('Write key {}'.format(key))

    #for key,item in outputs_metadict.items():
    #    f=key.replace(' ','_')
    #    item.to_pickle(f+'_meta.pkl')
    #    print('Write meta key {}'.format(key))

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--url', action='store', dest='url', default=None, type=str,
                        help='ASGS url to fetch ADCIRC data')
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    parser.add_argument('--data_source', action='store', dest='data_source', default='ASGS', type=str,
                        help='choose supported data source (case independant) eg ASGS')
    parser.add_argument('--data_product', action='store', dest='data_product', default='water_level', type=str,
                        help='choose supported data product eg water_level')
    parser.add_argument('--return_sample_min', action='store', dest='return_sample_min', default=60, type=int,
                        help='return_sample_min is the time stepping in the final data objects. (mins)')
    parser.add_argument('--ndays', default=-4, action='store', dest='ndays',help='Day lag (usually < 0)', type=int)
    parser.add_argument('--config_name', action='store', dest='config_name', default=None,
                        help='String: yml config which contains URL structural information')
    parser.add_argument('--contrails_auth', action='store', dest='contrails_auth', default=None, type=str,
                        help='Choose a non-default contrails auth config_name')
    parser.add_argument('--instance_name', action='store', dest='instance_name', default=None,
                        help='String: instance name')
    parser.add_argument('--fort63_style', action='store_true', default=False,
                        help='Boolean: Will inform Harvester to use fort.63.methods to get station nodesids')
    parser.add_argument('--hurricane_source', action='store',dest='hurricane_source', default=None,
                        help='str: Only needed for Hurricanes AND if using YAML to specify urls')
    parser.add_argument('--hurricane_year', action='store',dest='hurricane_year', default=None,
                        help='str: Only needed for Hurricanes AND if using YAML to specify urls')
    parser.add_argument('--gridname', action='store',dest='gridname', default='hsofs',
                        help='str: Select appropriate gridname Default is hsofs')
    parser.add_argument('--ensemble', action='store',dest='ensemble', default='nowcast',
                        help='str: Select appropriate ensemble Default is nowcast')
    parser.add_argument('--map_file', action='store',dest='map_file', default=None,
                        help='str: Select appropriate map_file yml for grid lookup')
    parser.add_argument('--main_yamlname', action='store',dest='main_yamlname', default=None,
                        help='str: Select appropriate main_yamlname')
    parser.add_argument('--finalDIR', action='store', dest='finalDIR', default=None,
                        help='String: Custom location for the output dicts, PNGs and logs')
    args = parser.parse_args()
    sys.exit(main(args))