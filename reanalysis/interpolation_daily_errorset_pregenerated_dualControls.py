#!/usr/bin/env python

# This code simply interpolates a single input data set

# The interpolator has been changed to LinearRBF. This requires an additional call to build a polygon boundary

# This code handles havong a separate water control line around PR

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

def station_list_interpolation(df, station_coords, station_list):
    """
    """
    model = interpolate_scaled_offset_field.interpolation_model_fit(df, interpolation_type='LinearRBF')
    df_interpolated = interpolate_scaled_offset_field.interpolation_model_transform(station_coords, model=model, input_grid_type='points',pathpoly=None)
    df_interpolated.index=station_list
    return df_interpolated

def df_adcirc_store(df_interpolated_ADCIRC_GRID, outputdir, iometadata=None):
    """
    """
    iosubdir='interpolated'
    gridfile = io_utilities.write_ADCIRC_formatted_gridfield_to_Disk(df_interpolated_ADCIRC_GRID, value_name='VAL', rootdir=outputdir,subdir='interpolated',fileroot='ADCIRC_interpolated_wl',iometadata=iometadata)
    utilities.log.info('Wrote ADCIRC offset field to {}'.format(gridfile))
    adcirc_interpolated_pkl = io_utilities.write_pickle(df_interpolated_ADCIRC_GRID, rootdir=outputdir,subdir=iosubdir,fileroot='interpolated_wl',iometadata=iometadata)
    utilities.log.info('Wrote ADCIRC offset field PKL {}'.format(adcirc_interpolated_pkl))
    print('Finished')

def model_store(model, outputdir, iometadata=None):
    """
    """    # Write out the model for posterity
    iosubdir='models'
    newfilename = io_utilities.get_full_filename_with_subdirectory_prepended(outputdir, iosubdir, 'interpolate_linear_model'+iometadata+'.h5')
    try:
        joblib.dump(model, newfilename)
        status = True
        utilities.log.info('Saved model file '+str(newfilename))
    except Exception as ex:
        utilities.log.error(f'Could not dump model file to disk {newfilename}. Ex {ex}')

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
## fetch the control point filenames
##
    file_land_controls = grid_to_station_maps.find_land_control_points_from_map(gridname=gridname, mapfile=args.map_file) 
    file_water_controls = grid_to_station_maps.find_water_control_points_from_map(gridname=gridname, mapfile=args.map_file)
## Find the secondary water controls for the island of PR
    pr_file_water_controls = grid_to_station_maps.find_secondary_water_control_points_from_map(gridname=gridname, mapfile=args.map_file)
    utilities.log.info(f'Fetched the clamp {file_land_controls}, and control node {file_water_controls} and secondary {pr_file_water_controls} files for grid {gridname}')
    print(f' secondary file {pr_file_water_controls}')

# Get the control points
    if file_water_controls is not None:
        df_land_controls = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=file_land_controls)
    else:
        df_land_controls=None
    if file_water_controls is not None:
        df_water_controls = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=file_water_controls)
    else:
        df_water_controls=None

    if pr_file_water_controls is not None:
        df_water_controls_pr = interpolate_scaled_offset_field.get_lon_lat_control_values(fname=pr_file_water_controls)
    else:
        df_water_controls_pr=None

    utilities.log.info(f'Control points {df_water_controls.shape}')
    utilities.log.info(f'Control points data {df_water_controls}')
    utilities.log.info(f'Secondary Control points data {df_water_controls_pr}')

# Get the station estimation data
    df_slr_stations = interpolate_scaled_offset_field.get_station_values(fname=args.slr_file_errors, header_data='VAL')
    df_mid_stations = interpolate_scaled_offset_field.get_station_values(fname=args.mid_file_errors, header_data='VAL')
    df_vlf_stations = interpolate_scaled_offset_field.get_station_values(fname=args.vlf_file_errors, header_data='VAL')

##
## Fetch the ADCIRC triangular grid file
##
#   griddir='/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/reanalysis/ec95d.v1/YEARLY-2019'
#   adc_coord_file = f'{griddir}/adc_coord.json'
    # Use this for HSOFS
    #adc_coords = io_utilities.read_json_file(f'{topdir}/adc_coord.json')
    #adc_grid_style = 'points' # Needed to guide the final interpolation
    #print(f' ADC_COORDS {adc_coords}')

##Fetch the EC95D coordinates

    griddir='/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/reanalysis/ec95d.v1/YEARLY-2019'
    adc_coord_file = f'{griddir}/adc_coord.json'
    adc_coords = io_utilities.read_json_file(adc_coord_file)
    adc_grid_style = 'points' # Needed to guide the final interpolation
    print(f' EC95D ADC_COORDS {adc_coords}')

##
## Fetch station_coords for which ALL data sets will be interpolated PRIOR to SUMMING
## As these are exact interpolators this is prob not necc for the MID/VF daa sets
##
    station_x = df_vlf_stations['LON']
    station_y = df_vlf_stations['LAT']
    station_list = df_vlf_stations.index
    station_coords = {'LON':station_x[:].tolist(),'LAT':station_y[:].tolist()}

##
## Prepare input data (SLR) to have the same input stations 
## No need to have CLAMPS for theaes model-specific interpolations
##
    df_slr_stations_interpolated=station_list_interpolation(df_slr_stations, station_coords, station_list)
    df_mid_stations_interpolated=station_list_interpolation(df_mid_stations, station_coords, station_list)
    df_vlf_stations_interpolated=station_list_interpolation(df_vlf_stations, station_coords, station_list)

    #model = interpolate_scaled_offset_field.interpolation_model_fit(df_slr_stations, interpolation_type='LinearRBF')
    #df_slr_stations_interpolated = interpolate_scaled_offset_field.interpolation_model_transform(station_coords, model=model, input_grid_type='points',pathpoly=None)
    #df_slr_stations_interpolated.index=df_vlf_stations.index

    #print(f' SLR input {df_slr_stations}')
    #print(f' SLR fit {df_slr_stations_interpolated}')

##
## Build the SUM data sets as mid_vlf-slr
##
    df_sum_stations_interpolated=df_vlf_stations_interpolated # MID and VLF always have the FINAL set of stations to keep
    df_sum_stations_interpolated['VAL']=df_mid_stations_interpolated['VAL']+df_vlf_stations_interpolated['VAL']-df_slr_stations_interpolated['VAL']
    print(f' SUMMED {df_sum_stations_interpolated}')
    print(f' SLR {df_slr_stations_interpolated}')

##
## Define new inteprolator function
##
    def df_model_surface_interpolation(df_stations, adc_coords, station_list, dataname=None, df_land_controls=None, df_water_controls=None, df_secondary_water_controls=None, cv_testing=None):
        """
        Dataname= (SLR,MID, or VLF)
        """
        in_interpolator='LinearRBF'
        set_of_dfs=list()
        set_of_poly_clamp_dfs=list()
        set_of_all_clamp_dfs=list()
    
        set_of_dfs.append(df_stations.rename(columns = {dataname:'VAL'})) # Always need this

        if df_water_controls is not None:
            set_of_dfs.append(df_water_controls)
            set_of_all_clamp_dfs.append(df_water_controls)
            set_of_poly_clamp_dfs.append(df_water_controls)
        if df_land_controls is not None:
            set_of_dfs.append(df_land_controls)
            set_of_all_clamp_dfs.append(df_land_controls)
            set_of_poly_clamp_dfs.append(df_land_controls)
        if df_secondary_water_controls is not None:
            set_of_dfs.append(df_secondary_water_controls)
            # Not used here set_of__clamp_dfs.append(df_secondary_water_controls)

        utilities.log.info('construct_interpolation_model: Number of dfs to combine for interpolation is {}'.format(len(set_of_dfs)))
        utilities.log.info('construct_interpolation_model: Number of dfs to combine for clamps polygon path building is {}'.format(len(set_of_all_clamp_dfs)))
    
        ## Use ALL data and clamps for interpoating and plotting
        #df_all_ClampControl = pd.concat(set_of_clamp_dfs,axis=0)
        #utilities.log.info(f'ClampControl {df_all_ClampControl}')

        # Buid separate polygons
        # This can only be used with LinearRBF. More testing is required
        df_poly_ClampControl=pd.concat(set_of_poly_clamp_dfs,axis=0)
        df_poly_pr_ClampControl=df_secondary_water_controls

        print('secondary {df_poly_pr_ClampControl}')
        pathpoly = interpolate_scaled_offset_field.build_polygon_path(df_poly_ClampControl) 
        pathpoly_pr = interpolate_scaled_offset_field.build_polygon_path(df_poly_pr_ClampControl)
        df_combined=interpolate_scaled_offset_field.combine_datasets_for_interpolation(set_of_dfs)

        model = interpolate_scaled_offset_field.interpolation_model_fit(df_combined, interpolation_type=in_interpolator)

        # Build new extrapolated surface onto the ADCIRC grid
        df_interpolated_ADCIRC_GRID = interpolate_scaled_offset_field.interpolation_model_transform_dualpoly(adc_coords, model=model, input_grid_type='points',pathpoly=pathpoly, pathpoly2=pathpoly_pr)

        #if cv_testing:
        #    full_scores, best_scores = interpolate_scaled_offset_field.test_interpolation_fit(df_source=df_stations[['LON','LAT',dataname]], df_land_controls=df_land_controls, df_water_controls=df_water_controls, cv_splits=5)
        #    print(best_scores)
        #    print(full_scores)
        #    combined_scoring={'best_scores':best_scores, 'full_scores':full_scores, 'Interpolator':in_interpolator}
        #    utilities.log.info('Performed a CV testing')

        # For subsequent data viewers change the VAL header back to the real data name
        return model, df_interpolated_ADCIRC_GRID

##
## Process the four clamps Interpolation surfaces
##

# NEED to send in the secondary water controls

    model_slr,df_interpolated_ADCIRC_GRID_SLR=df_model_surface_interpolation(df_slr_stations_interpolated, adc_coords, station_list, dataname='SLR', df_land_controls=df_land_controls, df_water_controls=df_water_controls, 
				df_secondary_water_controls=df_water_controls_pr, cv_testing=None)

    model_mid,df_interpolated_ADCIRC_GRID_MID=df_model_surface_interpolation(df_mid_stations_interpolated, adc_coords, station_list, dataname='SLR', df_land_controls=df_land_controls, df_water_controls=df_water_controls, 
				df_secondary_water_controls=df_water_controls_pr, cv_testing=None)

    model_vlf,df_interpolated_ADCIRC_GRID_VLF=df_model_surface_interpolation(df_vlf_stations_interpolated, adc_coords, station_list, dataname='SLR', df_land_controls=df_land_controls, df_water_controls=df_water_controls, 
				df_secondary_water_controls=df_water_controls_pr, cv_testing=None)

    model_sum,df_interpolated_ADCIRC_GRID_SUM=df_model_surface_interpolation(df_sum_stations_interpolated, adc_coords, station_list, dataname='SLR', df_land_controls=df_land_controls, df_water_controls=df_water_controls, 
				df_secondary_water_controls=df_water_controls_pr, cv_testing=None)
    print('After interpolations')
    print(f' shit1 {df_slr_stations_interpolated}')
    print(f' shit2 {df_sum_stations_interpolated}')

##
## Write out datafiles for the chgosen models
##

    df_adcirc_store(df_interpolated_ADCIRC_GRID_SLR, outputdir, iometadata=f'{iometadata}_SLR')
    df_adcirc_store(df_interpolated_ADCIRC_GRID_MID, outputdir, iometadata=f'{iometadata}_MID')
    df_adcirc_store(df_interpolated_ADCIRC_GRID_VLF, outputdir, iometadata=f'{iometadata}_VLF')
    df_adcirc_store(df_interpolated_ADCIRC_GRID_SUM, outputdir, iometadata=f'{iometadata}_SUM')
    model_store(model_slr, outputdir, iometadata=f'{iometadata}_SLR')
    model_store(model_mid, outputdir, iometadata=f'{iometadata}_MID')
    model_store(model_vlf, outputdir, iometadata=f'{iometadata}_VLF')
    model_store(model_sum, outputdir, iometadata=f'{iometadata}_SUM')
    print('HANDLED the IO')

    # If scoring CV, then write out that data as well.
    #iosubdir='interpolated'
    #if cv_testing:
    #     score_filename=io_utilities.write_dict_to_json(combined_scoring, rootdir=outputdir,subdir=iosubdir,fileroot='reanalysis_CV_scoring',iometadata=iometadata)
    #     utilities.log.info('Wrote CV scoring data to {}'.format(score_filename))

    # For posterity also dump out the KNN values of the processed control nodes ( if provided)

    #if df_land_controls is not None:
    #    iosubdir='interpolated'
    #    file_land_controls = io_utilities.write_json(df_land_controls, rootdir=outputdir,subdir=iosubdir,fileroot='controlNodeSummary',iometadata=iometadata) 
    #    utilities.log.info('Processed land control nodes: Final values stored to {}'.format(file_land_controls))

####################################################################################################################

## Yes we want this as well
##
## Optional. Apply the model to a 500x400 grid and plot, the extrapolated surface, stations, clamps
##
## Just use the SUM for now

    # Test using the generic grid and plot to see the generated offset surface. Use the same model as previously generated
    #adc_plot_grid = interpolate_scaled_offset_field.generic_grid()
    #df_plot_transformed = interpolate_scaled_offset_field.interpolation_model_transform_dualpoly(adc_plot_grid, model=model_sum, input_grid_type='grid', pathpoly=pathpoly, pathpoly2=pathpoly_pr)
    #iosubdir='images'
    #annotate = f'{gridname.upper()}_{iometadata}'
    #newfilename = io_utilities.get_full_filename_with_subdirectory_prepended(outputdir, iosubdir, 'interpolated_surface_plot_'+iometadata+'.png')
    ##df_sum_stations_interpolated
    #df_all_water_controls = pd,concat([df_water_controls,df_secondary_water_controls])

    #plot_interpolation_errorset.save_plot_model( adc_plot_grid=adc_plot_grid, df_surface=df_plot_transformed, df_stations=df_sum_stations, df_land_control=df_land_controls, df_all_water_control=df_all_water_controls, filename=newfilename, plot_now=False, annotation=annotate)
    #utilities.log.info('Saved IMAGE file to {}'.format(newfilename))
    #utilities.log.info('Finished with interpolation')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--main_yamlname', action='store',dest='main_yamlname', default=None,
                        help='str: Select appropriate main_yamlname')
    parser.add_argument('--slr_file_errors', action='store', dest='slr_file_errors', default=None,
                        help='CSV that contains the compute_error The current station list')
    parser.add_argument('--mid_file_errors', action='store', dest='mid_file_errors', default=None,
                        help='CSV that contains the compute_error The current station list')
    parser.add_argument('--vlf_file_errors', action='store', dest='vlf_file_errors', default=None,
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
