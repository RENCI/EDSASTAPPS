#!/usr/bin/env python
#

import os,sys
import numpy as np
import pandas as pd
import datetime as dt

import harvester.fetch_adcirc_data as fetch_adcirc_data
import harvester.generate_urls_from_times as genurls # generate_urls_from_times
import harvester.get_adcirc_stations as get_adcirc_stations
import harvester.get_observations_stations as get_obs_stations
import processing.compute_error_field as compute_error_field
import processing.interpolate_scaled_offset_field as interpolate_scaled_offset_field
import gridmap.grid_to_station_maps as grid_to_station_maps
import adda_visualization_plots as adda_visualization_plots

#from harvester.fetch_adcirc_data import fetch_adcirc_data
#from harvester.generate_urls_from_times import genurls # generate_urls_from_times
#from harvester.get_adcirc_stations import get_adcirc_stations
#from harvester.get_observations_stations import get_obs_stations
#from processing.compute_error_field import compute_error_field
#from processing.interpolate_scaled_offset_field import interpolate_scaled_offset_field
#from grid_to_station_maps import grid_to_station_maps
#from adda_visualization_plots import adda_visualization_plots

from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities
from argparse import ArgumentParser

import joblib

##
## Run the canonical ADDA approach
##

# user must coordinate the proper AST data_source (NOAA< NDBC, etc) to the chosen station list
 
def main(args):
    """
    A simple main method to demonstrate the use of this code
    It assumes the existance of a proper main.yml to get IO information
    It assumes the existance of a proper url_framework.yml which (optionally) can be used to create URLs
    It assumes the existance of a proper local_instance.yml which (optionally_) may be used to read customized annual Instance values
    """

    config = utilities.init_logging(subdir=None, config_file='./config/main.yml')

    # Basic checks
    if args.config_name is None:
        config_name =os.path.join(os.path.dirname(__file__), './secrets', 'url_framework.yml')
    if args.instance_name is None:
        utilities.log.error('Instance name is missing but is mandatory')
        sys.exit(1)

    hurricane_source=args.hurricane_source
    hurricane_year=args.hurricane_year

    fort63_style=args.fort63_style

    # Set up IO env
    utilities.log.info("Product Level Working in {}.".format(os.getcwd()))

    # Set up the times information No need to worry about values for hh:mm:ssZ Subsequent resampling cleans that up
    if args.timeout is None: # Set to a default of now()
        tnow = dt.datetime.now()
        stoptime = tnow.strftime('%Y-%m-%d %H:%M:%S')
    else:
        stoptime=args.timeout
    print('Stoptime and ndays {}. {}'.format(stoptime,args.ndays))

##
## fetch the asgs/adcirc station data
##
    #map_file=os.path.join(os.path.dirname(__file__), './supporting_data', 'grid_to_stationfile_maps.yml')
    #Note specifying the map_file REQUIRES the listed file to have a fullpathname
    station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=args.gridname, mapfile=args.map_file, datatype='NOAA_STATIONS')
    file_land_controls = grid_to_station_maps.find_land_control_points_from_map(gridname=args.gridname, mapfile=args.map_file)   # ,mapfile=map_file)
    file_water_controls = grid_to_station_maps.find_water_control_points_from_map(gridname=args.gridname, mapfile=args.map_file)

    # Specify file containing the list of stations
    if fort63_style:
        if not fort63_compliant:
            utilities.log.error('Selected grid list is not fort63 compliant')
            sys.exit(1)

    # Generate a list of URLs consistent with the desired criteria
    print('Demonstrate URL generation spanning multiple instance values (hsofs). YAML, timeout, offset')

    #hurricane_source=None
    #hurricane_year=None
    rpl = genurls.generate_urls_from_times(timeout=stoptime, ndays=args.ndays, grid_name=args.gridname, instance_name=args.instance_name, config_name=config_name, hurricane_yaml_year=hurricane_year,hurricane_yaml_source=hurricane_source)
    urls = rpl.build_url_list_from_yaml_and_offset(ensemble=args.ensemble)

    print(station_file)

    rpl = get_adcirc_stations.get_adcirc_stations(source='ASGS', product=args.data_product,
                station_list_file=station_file, 
                knockout_file=None, fort63_style=fort63_style )

    # Convert URLs to desired fort type
    if fort63_style:
        urls=get_adcirc_stations.convert_urls_to_63style(urls)
    else:
        urls=get_adcirc_stations.convert_urls_to_61style(urls)
    print(urls)

    # Fetch best resolution and no resampling
    data_adc,meta_adc=rpl.fetch_station_product(urls, return_sample_min=args.return_sample_min, fort63_style=fort63_style  )

    # Revert Harvester filling of nans to -99999 back to nans
    data_adc.replace('-99999',np.nan,inplace=True)
    meta_adc.replace('-99999',np.nan,inplace=True)

    # Grab the grid coordinates for the url 
    urls_63 = get_adcirc_stations.convert_urls_to_63style(urls)
    print(urls_63)
    adc_coords = get_adcirc_stations.extract_adcirc_grid_coords( urls_63 )

    #lons = adc_coords['LON']
    #lats = adc_coords['LAT']
    #adc_coords = {'lon':lons, 'lat':lats} # This created for backward compatibility

    print(f'Grid name {rpl.gridname}')
    print(f'Instance name {rpl.instance}')

    # Get a last piece of metadata for firsr url in iterable grabs either the time (%Y%m%s%H) or hurricane advisory (int)
    # Grab the adcirc time ranges for calling the observations code
    
    obs_starttime = dt.datetime.strftime( min(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')
    obs_endtime = dt.datetime.strftime( max(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')

    # Cnstruct iometadata to update all filename - want diff format
    io_start = dt.datetime.strftime( min(data_adc.index.tolist()), '%Y%m%d%H%M')
    io_end = dt.datetime.strftime( max(data_adc.index.tolist()), '%Y%m%d%H%M')

    iometadata = io_start+'_'+io_end 
    rootdir=io_utilities.construct_base_rootdir(config['DEFAULT']['RDIR'], base_dir_extra='ADDA'+'_'+iometadata)
    
    # Write the data to disk in a way that mimics ADDA

    iosubdir='adcpkl'
    # Write selected in Pickle data 
    metapkl = io_utilities.write_pickle(meta_adc,rootdir=rootdir,subdir=iosubdir,fileroot='adc_wl_metadata',iometadata=iometadata)
    detailedpkl = io_utilities.write_pickle(data_adc, rootdir=rootdir,subdir=iosubdir,fileroot='adc_wl_detailed',iometadata=iometadata)

    # Grab the adcirc time ranges for calling the observations code
    obs_starttime = dt.datetime.strftime( min(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')
    obs_endtime = dt.datetime.strftime( max(data_adc.index.tolist()), '%Y-%m-%d %H:%M:%S')

    # Optionally Write selected in JSON format Convert and write selected JSON data
    #data_adc.index = data_adc.index.strftime('%Y-%m-%d %H:%M:%S')
    data_adc_4json = data_adc.copy()
    data_adc_4json.index = data_adc_4json.index.strftime('%Y-%m-%d %H:%M:%S')
    metajson = io_utilities.write_json(meta_adc,rootdir=rootdir,subdir=iosubdir,fileroot='adc_wl_metadata',iometadata=iometadata)
    detailedjson = io_utilities.write_json(data_adc_4json, rootdir=rootdir,subdir=iosubdir,fileroot='adc_wl_detailed',iometadata=iometadata)
   
    print('Write out coordinates')
    #ADCfilecoords=utilities.write_csv(df_adc_coords, rootdir=rootdir,subdir=iosubdir,fileroot='adc_coord',iometadata=iometadata)
    #utilities.log.info('Wrote grid coords to {}'.format(ADCfilecoords))
    print('Write out coordinates dict in JSON format')
    ADCfilecoordsJson = io_utilities.write_dict_to_json(adc_coords, rootdir=rootdir,subdir=iosubdir,fileroot='adc_coord',iometadata=iometadata)
    #ADCfilecoordsJson = rootdir+'/adc_coord'+iometadata+'.json'
    io_utilities.write_json_file(adc_coords, ADCfilecoordsJson)
    utilities.log.info('Wrote grid coords to {}'.format(ADCfilecoordsJson))

##
## set up and fetch the NOAA observations - using a time range specified by the available adcirc data set
##

    # Could try to use contrails, too.
    #contrails_stations=fname=os.path.join(os.path.dirname(__file__),'./supporting_data','contrails_stations_coastal.csv')
    #contrails_stations=fname=os.path.join(os.path.dirname(__file__),'./supporting_data','contrails_stations_rivers.csv')

    obs = get_obs_stations.get_obs_stations(source='NOAA', product='water_level',
                contrails_yamlname='None',
                knockout_file=None, station_list_file=station_file)
    # Get data at highest resolution
    data_obs,meta_obs=obs.fetch_station_product((obs_starttime,obs_endtime), return_sample_min=0, interval='None' )
    data_obs.replace('-99999',np.nan,inplace=True)
    meta_obs.replace('-99999',np.nan,inplace=True)
    # Remove stations with too many nans ( Note Harvester would have previously removed stations that are ALL NANS)
    data_thresholded = obs.remove_missingness_stations(data_obs, max_nan_percentage_cutoff=10)  # (Maximum allowable nans %)
    # meta_thresholded = meta_obs.loc[data_thresholded.columns.tolist()]
    meta = set(data_thresholded.columns.tolist()).intersection(meta_obs.index.to_list())
    meta_thresholded = meta_obs.loc[meta]
    # Apply a moving average (smooth) the data performed the required resampling to the desired rate followed by interpolating
    data_obs_smoothed = obs.fetch_smoothed_station_product(data_thresholded, return_sample_min=60, window=11)

    # Write the data to disk in a way that mimics ADDA
    iosubdir='obspkl'
    # Write selected in Pickle data 
    metapkl = io_utilities.write_pickle(meta_thresholded,rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_metadata',iometadata=iometadata)
    detailedpkl = io_utilities.write_pickle(data_thresholded, rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_detailed',iometadata=iometadata)
    smoothpkl = io_utilities.write_pickle(data_obs_smoothed, rootdir=rootdir,subdir=iosubdir,fileroot='obs_wl_smoothed',iometadata=iometadata)
    print('Finished')

##
## Perform the error computations
##

    #bound_lo='2022-02-17 06:00:00'
    #bound_hi='2022-02-19 00:00:00'
    print(data_obs_smoothed.index)
    comp_err = compute_error_field.compute_error_field(data_obs_smoothed, data_adc, meta_obs) # All default params
    comp_err._intersection_stations()
    comp_err._intersection_times()
    comp_err._tidal_transform_data()
    comp_err._apply_time_bounds((obs_starttime,obs_endtime)) # redundant but here for illustration
    comp_err._compute_and_average_errors()

    # Set up IO env
    utilities.log.info("Product Level Working in {}.".format(os.getcwd()))
    main_yamlname=os.path.join(os.path.dirname(__file__), './config', 'main.yml')
    config = utilities.load_config(main_yamlname)
    #rootdir=utilities.fetchBasedir(config['DEFAULT']['RDIR'], basedirExtra='ADDA'+iometadata)

    iosubdir='errorfield'

    # Write selected Pickle data 
    tideTimeErrors = io_utilities.write_pickle(comp_err.diff,rootdir=rootdir,subdir=iosubdir,fileroot='tideTimeErrors',iometadata=iometadata)
    utilities.log.info('Wrote out Pickle {}'.format(tideTimeErrors))

    # Write selected in JSON format. Basically the Merged_dict data multi-index set
    data_dict = compute_error_field.combine_data_to_dict(comp_err.adc,comp_err.obs,comp_err.diff, product='WL')
    dict_json = io_utilities.write_dict_to_json(data_dict, rootdir=rootdir,subdir=iosubdir,fileroot='adc_obs_error_merged',iometadata=iometadata)
    utilities.log.info('Wrote out JSON {}'.format(dict_json))

    # Write the CSV files which carry averages
    station_summary_aves=io_utilities.write_csv(comp_err.df_final, rootdir=rootdir,subdir=iosubdir,fileroot='stationSummaryAves',iometadata=iometadata)
    station_period_aves=io_utilities.write_csv(comp_err.df_cycle_aves, rootdir=rootdir,subdir=iosubdir,fileroot='stationPeriodAves',iometadata=iometadata)
    utilities.log.info('Wrote out CSV cyclic averages {}'.format(station_period_aves))
    utilities.log.info('Wrote out CSV summary averages {}'.format(station_summary_aves))

##
## Perform the interpolation
##
    # map_file=os.path.join(os.path.dirname(__file__), './supporting_data', 'grid_to_stationfile_maps.yml')
    file_stations = station_summary_aves # Better to reread using "get_station_values"
    #file_land_controls = grid_to_station_maps.read_specified_filename_from_map(gridname=args.gridname, datatype='LANDCONTROL')   #  mapfile=map_file
    #file_water_controls = grid_to_station_maps.read_specified_filename_from_map(gridname=args.gridname, datatype='WATERCONTROL') #  mapfile=map_file
    #
    df_stations = interpolate_scaled_offset_field.get_station_values(fname=file_stations)
    df_land_controls = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=file_land_controls)
    df_water_controls = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=file_water_controls)

    set_of_dfs=list()
    set_of_dfs.append(df_stations) # Always need this
    if df_land_controls is not None:
        # Must remover NaNs
        df_new_land_controls = interpolate_scaled_offset_field.knn_fit_control_points(df_stations, df_land_controls, nearest_neighbors=3)
        set_of_dfs.append(df_new_land_controls)
    if df_water_controls is not None:
        set_of_dfs.append(df_water_controls)
    utilities.log.info('construct_interpolation_model: Number of dfs to combine for interpolation is {}'.format(len(set_of_dfs)))
    df_combined=interpolate_scaled_offset_field.combine_datasets_for_interpolation(set_of_dfs)
    model = interpolate_scaled_offset_field.interpolation_model_fit(df_combined, fill_value=0.0, interpolation_type='LinearNDInterpolator')

    # Build new grid
    df_extrapolated_ADCIRC_GRID = interpolate_scaled_offset_field.interpolation_model_transform(adc_coords, model=model, input_grid_type='points')

    # do A TEST FIT
    print('TEST FIT')
    full_scores, best_scores = interpolate_scaled_offset_field.test_interpolation_fit(df_source=df_stations, df_land_controls=df_land_controls, df_water_controls=df_water_controls, cv_splits=5, nearest_neighbors=3)
    print('ADDA')
    print(best_scores)
    print(full_scores)
##
## Write out datafiles 
##
    gridfile = io_utilities.write_ADCIRC_formatted_gridfield_to_Disk(df_extrapolated_ADCIRC_GRID, value_name='VAL', rootdir=rootdir,subdir='interpolated',fileroot='ADCIRC_interpolated_wl',iometadata=iometadata)
    utilities.log.info('Wrote ADCIRC offset field to {}'.format(gridfile))
    adcirc_extrapolated_pkl = io_utilities.write_pickle(df_extrapolated_ADCIRC_GRID, rootdir=rootdir,subdir=iosubdir,fileroot='interpolated_wl',iometadata=iometadata)
    utilities.log.info('Wrote ADCIRC offset field PKL {}'.format(adcirc_extrapolated_pkl))
    print('Finished')

    # Test using the generic grid and plot to see the generated offset surface. Use the same model as previously generated
    adc_plot_grid = interpolate_scaled_offset_field.generic_grid()
    df_plot_transformed = interpolate_scaled_offset_field.interpolation_model_transform(adc_plot_grid, model=model, input_grid_type='grid') 

    # Write out the model for posterity
    newfilename = io_utilities.get_full_filename_with_subdirectory_prepended(rootdir, iosubdir, 'interpolate_linear_model'+iometadata+'.h5')
    try:
        joblib.dump(model, newfilename)
        status = True
        utilities.log.info('Saved model file '+str(newfilename))
    except:
        utilities.log.error('Could not dump model file to disk '+ newfilename)

##
## Optional. Apply the model to a 400x300 grid and plot, the extrapolated surface, stations, clamps
##

# Choose to plot or save the plot file ?

    iosubdir='images'
    newfilename = io_utilities.get_full_filename_with_subdirectory_prepended(rootdir, iosubdir, 'extrapolated_surface_plot'+iometadata+'.png')
    adda_visualization_plots.save_plot_model( adc_plot_grid=adc_plot_grid, df_surface=df_plot_transformed, df_land_control=df_land_controls, df_water_control=df_water_controls, filename=newfilename, plot_now=False)
    utilities.log.info('Saved IMAGE file to {}'.format(newfilename))

##
## Dump temp files for use in testing
##
    
    # Write out temp files for testing visualization
    #df_plot_transformed.to_csv('/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_surface.csv')
    #df_new_land_controls.to_csv('/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_land_controls.csv')
    #df_water_controls.to_csv('/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2/test_data/df_water_controls.csv')
    #ADCJson = io_utilities.write_dict_to_json(adc_plot_grid, rootdir='/projects/sequence_analysis/vol1/prediction_work/ADCIRCSupportTools-v2',subdir='test_data',fileroot='adc_plot_grid',iometadata='')




# We need to support both specifying URLs by explicit urls and by specifying time ranges.

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--sources', action='store_true',
                        help='List currently supported data sources')
    #parser.add_argument('--data_source', action='store', dest='data_source', default='ASGS', type=str,
    #                    help='choose supported data source (case independant) eg ASGS')
    parser.add_argument('--data_product', action='store', dest='data_product', default='water_level', type=str,
                        help='choose supported data product eg water_level')
    parser.add_argument('--return_sample_min', action='store', dest='return_sample_min', default=60, type=int,
                        help='return_sample_min is the time stepping in the final data objects. (mins)')
    parser.add_argument('--ndays', default=-2, action='store', dest='ndays',help='Day lag (usually < 0)', type=int)
    parser.add_argument('--timeout', default=None, action='store', dest='timeout', help='YYYY-mm-dd HH:MM:SS. Latest day of analysis def to now()', type=str)
    parser.add_argument('--config_name', action='store', dest='config_name', default=None,
                        help='String: yml config which contains URL structural information')
    parser.add_argument('--instance_name', action='store', dest='instance_name', default=None,
                        help='String: instance name')
    parser.add_argument('--fort63_style', action='store_true',
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
                        help='str: Select appropriate map_file ym; for grid lookup')
    args = parser.parse_args()
    sys.exit(main(args))