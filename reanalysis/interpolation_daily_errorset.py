#!/usr/bin/env python

# This code simply interpolates a single input data set

# The interpolator has been changed to LinearRBF. This requires an additional call to build a polygon boundary

###########################################################################

## Added code to capture the set of scores for subsequenct analysis

import os,sys
import numpy as np
import pandas as pd
import time as tm
import datetime as dt
import processing.interpolate_scaled_offset_field as interpolate_scaled_offset_field
import gridmap.grid_to_station_maps as grid_to_station_maps
import plot_interpolation_errorset as plot_interpolation_errorset
import joblib

from utilities.utilities import utilities as utilities
import io_utilities.io_utilities as io_utilities

###########################################################################
    
def fetch_data_metadata(metaFile=None):
    metaDict=None
    try:
        #meta='/'.join([metaFile,'obs_water_level_metadata.json'])
        with open(meta, 'r') as fp:
            try:
                metaDict = json.load(fp)
            except OSError:
                utilities.log.error("Could not open/read file {}".format(meta))
                sys.exit()
    except TypeError:
        utilities.log.info('No directory specified. Set metadata filename to None')
    return metaDict

def main(args):
    #assumes the existance of a proper main.yml to get IO information
    if args.main_yamlname is None:
        main_yamlname=os.path.join(os.path.dirname(__file__), './config', 'main.yml')
    else:
        main_yamlname=args.main_yamlname
    config = utilities.init_logging(subdir=None, config_file=main_yamlname, log_file_metadata='interpolate')
    utilities.log.info('Selected main_yamlfile of {}'.format(main_yamlname))

    if not args.input_directory:
        utilities.log.error('Need input_directory on command line: --input_directory <input_directory>')
        return 1
    topdir = args.input_directory.strip()

    if not args.output_directory:
        utilities.log.error('Need output_directory on command line: --input_directory <input_directory>')
        return 1

    outputdir=io_utilities.construct_base_rootdir(args.output_directory.strip(), base_dir_extra='')

    cv_testing = args.cv_testing
    gridname = args.gridname
    iometadata = args.iometadata

##
## fetch the asgs/adcirc station filenames
##
    ##station_file, fort63_compliant = grid_to_station_maps.find_station_list_from_map(gridname=gridname, mapfile=args.map_file, datatype='NOAA_STATIONS')
    file_land_controls = grid_to_station_maps.find_land_control_points_from_map(gridname=gridname, mapfile=args.map_file) 
    file_water_controls = grid_to_station_maps.find_water_control_points_from_map(gridname=gridname, mapfile=args.map_file)
    utilities.log.info('Fetched the clamp {}, and control node {} files for grid {}'.format(file_land_controls,file_water_controls,gridname))

# Get the data
    df_stations = interpolate_scaled_offset_field.get_station_values(fname=args.daily_file_errors, header_data='fft')
    if file_land_controls is not None:
        df_land_controls = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=file_land_controls)
    else:
        df_land_controls=None
    if file_water_controls is not None:
        df_water_controls = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=file_water_controls)
    else:
        df_water_controls=None

    #utilities.log.info(f'Control points {df_water_controls.shape}')
    #utilities.log.info(f'Control points data {df_water_controls}')
    #utilities.log.info(f'Station input file is {args.daily_file_errors}')
    #utilities.log.info(f'Stations {df_stations}')

##
## Fetch the ADCIRC triangular grid file
##
    adc_coords = io_utilities.read_json_file(f'{topdir}/adc_coord.json')
    adc_grid_style = 'points' # Needed to guide the final interpolation

##
## Perform the Interpolation
##
    in_interpolator='LinearRBF'
    ##in_interpolator='NearestNDInterpolator'
    #in_nearest_neighbors=3
    set_of_dfs=list()
    set_of_clamp_dfs=list()
    set_of_dfs.append(df_stations) # Always need this
    if file_water_controls is not None:
        set_of_dfs.append(df_water_controls)
        set_of_clamp_dfs.append(df_water_controls)
    if file_land_controls is not None:
        utilities.log.info('No KNN fit the control points')
        set_of_dfs.append(df_land_controls)
        set_of_clamp_dfs.append(df_land_controls)
        utilities.log.info(f'Final Land control results {df_land_controls}')
    utilities.log.info('construct_interpolation_model: Number of dfs to combine for interpolation is {}'.format(len(set_of_dfs)))
    utilities.log.info('construct_interpolation_model: Number of dfs to combine for clamps polygon path building is {}'.format(len(set_of_clamp_dfs)))

    df_ClampControl = pd.concat(set_of_clamp_dfs,axis=0)
    utilities.log.info(f'ClampControl {df_ClampControl}')

    # This can only be used with LinearRBF. More testing is required
    pathpoly = interpolate_scaled_offset_field.build_polygon_path(df_ClampControl)

    df_combined=interpolate_scaled_offset_field.combine_datasets_for_interpolation(set_of_dfs)
    #utilities.log.info(f' Interpolation post KNN {df_combined.shape}')
    #utilities.log.info(f'After KNN: Stations {df_stations}')
    model = interpolate_scaled_offset_field.interpolation_model_fit(df_combined, interpolation_type=in_interpolator)

    # Apply model to the input data as a test. For non-exact interpolation methods this may be informative
    station_x = df_stations['LON']
    station_y = df_stations['LAT']
    station_coords = {'LON':station_x[:].tolist(),'LAT':station_y[:].tolist()}
    utilities.log.info(f'Interpolation. Number of input lon/lat pairs is {station_coords}')
    df_interpolated_stations = interpolate_scaled_offset_field.interpolation_model_transform(station_coords, model=model, input_grid_type='points',pathpoly=pathpoly)

    df_stations['interpolated']=df_interpolated_stations['VAL'].tolist()

    # Build new extrapolated surface onto the ADCIRC grid
    df_extrapolated_ADCIRC_GRID = interpolate_scaled_offset_field.interpolation_model_transform(adc_coords, model=model, input_grid_type='points',pathpoly=pathpoly)

    # do A TEST FIT and include other information for posterity. NOTE: under the covers, controls get merged with df_source. df_source must only enter
    # with the columns lON,LAT,VAL

    if cv_testing:
        full_scores, best_scores = interpolate_scaled_offset_field.test_interpolation_fit(df_source=df_stations[['LON','LAT','VAL']], df_land_controls=df_land_controls, df_water_controls=df_water_controls, cv_splits=5)
        print(best_scores)
        print(full_scores)
        combined_scoring={'best_scores':best_scores, 'full_scores':full_scores, 'Interpolator':in_interpolator}
        #combined_scoring={'best_scores':best_scores, 'full_scores':full_scores, 'KNN':in_nearest_neighbors, 'Interpolator':in_interpolator}
        utilities.log.info('Performed a CV testing')

    # For subsequent data viewers change the VAL header back to the real data name
    df_stations=df_stations.rename(columns = {'VAL':'fft'})

##
## Write out datafiles 
##
    iosubdir='interpolated'
    gridfile = io_utilities.write_ADCIRC_formatted_gridfield_to_Disk(df_extrapolated_ADCIRC_GRID, value_name='VAL', rootdir=outputdir,subdir='interpolated',fileroot='ADCIRC_interpolated_wl',iometadata=iometadata)
    utilities.log.info('Wrote ADCIRC offset field to {}'.format(gridfile))
    adcirc_extrapolated_pkl = io_utilities.write_pickle(df_extrapolated_ADCIRC_GRID, rootdir=outputdir,subdir=iosubdir,fileroot='interpolated_wl',iometadata=iometadata)
    utilities.log.info('Wrote ADCIRC offset field PKL {}'.format(adcirc_extrapolated_pkl))
    print('Finished')

    # Test using the generic grid and plot to see the generated offset surface. Use the same model as previously generated
    adc_plot_grid = interpolate_scaled_offset_field.generic_grid()
    df_plot_transformed = interpolate_scaled_offset_field.interpolation_model_transform(adc_plot_grid, model=model, input_grid_type='grid', pathpoly=pathpoly)

    # Write out the model for posterity
    iosubdir='models'
    newfilename = io_utilities.get_full_filename_with_subdirectory_prepended(outputdir, iosubdir, 'interpolate_linear_model'+iometadata+'.h5')
    try:
        joblib.dump(model, newfilename)
        status = True
        utilities.log.info('Saved model file '+str(newfilename))
    except:
        utilities.log.error('Could not dump model file to disk '+ newfilename)

    # Write out the new df_stations that includes the interpolate column
    # must build a new name from the old convert from stationSummaryAves_00-366_2000123100.csv to interpolate_stationSummary_00-366_2000123100.csv

    utilities.log.info(f'Write out interpolation data {df_stations.shape}')
    iosubdir='interpolated'
    new_sum_avefile = io_utilities.write_csv(df_stations, rootdir=outputdir,subdir=iosubdir,fileroot='interpolated_stationSummaryAves',iometadata=iometadata)
    utilities.log.info('Wrote Offset field CSV that includes interpolated values {}'.format(new_sum_avefile))

    # If scoring CV, then write out that data as well.
    iosubdir='interpolated'
    if cv_testing:
         score_filename=io_utilities.write_dict_to_json(combined_scoring, rootdir=outputdir,subdir=iosubdir,fileroot='reanalysis_CV_scoring',iometadata=iometadata)
         utilities.log.info('Wrote CV scoring data to {}'.format(score_filename))

    # For posterity also dump out the KNN values of the processed control nodes ( if provided)

    if df_land_controls is not None:
        iosubdir='interpolated'
        file_land_controls = io_utilities.write_json(df_land_controls, rootdir=outputdir,subdir=iosubdir,fileroot='controlNodeSummary',iometadata=iometadata) 
        utilities.log.info('Processed land control nodes: Final values stored to {}'.format(file_land_controls))

##
## Optional. Apply the model to a 500x400 grid and plot, the extrapolated surface, stations, clamps
##
    iosubdir='images'
    annotate = f'{gridname.upper()}_{iometadata}'
    newfilename = io_utilities.get_full_filename_with_subdirectory_prepended(outputdir, iosubdir, 'extrapolated_surface_plot_'+iometadata+'.png')
    plot_interpolation_errorset.save_plot_model( adc_plot_grid=adc_plot_grid, df_surface=df_plot_transformed, df_stations=df_stations.rename(columns={"fft": "VAL"}), df_land_control=df_land_controls, df_water_control=df_water_controls, filename=newfilename, plot_now=False, annotation=annotate)

    # Can save df_extrapolated_ADCIRC_GRID as a pkl with headers LON,LAT,VAL
    # This is a large memory operation but here is how the ADCIRC data would be be plotted
    # plot_interpolation_errorset.save_plot_model( adc_plot_grid=df_extrapolated_ADCIRC_GRID, df_surface=df_plot_transformed, df_stations=df_stations.rename(columns={"fft": "VAL"}), df_land_control=df_new_land_controls, df_water_control=df_water_controls, filename=newfilename, plot_now=False, annotation=annotate)

    utilities.log.info('Saved IMAGE file to {}'.format(newfilename))
    utilities.log.info('Finished with interpolation')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--main_yamlname', action='store',dest='main_yamlname', default=None,
                        help='str: Select appropriate main_yamlname')
    parser.add_argument('--daily_file_errors', action='store', dest='daily_file_errors', default=None,
                        help='CSV that contains the compute_error The current station list')
    parser.add_argument('--input_directory', action='store', dest='input_directory', default=None,
                        help='directory for yearly data')
    parser.add_argument('--output_directory', action='store', dest='output_directory', default=None,
                        help='directory for yearly data')
    parser.add_argument('--gridname', action='store',dest='gridname', default='hsofs',
                        help='str: Select appropriate gridname Default is hsofs')
    parser.add_argument('--cv_testing', action='store_true', dest='cv_testing',
                        help='Boolean: Invoke a CV procedure prior to fitting kriging model')
    parser.add_argument('--map_file', action='store',dest='map_file', default=None,
                        help='str: Select appropriate map_file ym; for grid lookup')
    parser.add_argument('--iometadata', action='store',dest='iometadata', default=None,
                        help='str: Added with all output files')

    args = parser.parse_args()
    sys.exit(main(args))
