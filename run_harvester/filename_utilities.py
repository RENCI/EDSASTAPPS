##############
## Methods to estimate what the output list of data and metadata filenames are supposed to be. We needs this
## because it is required as an argument to the Pegasus Jobs object specified before executing the a actual script
##############

import glob, os
import yaml
import json
import shutil
from pathlib import Path
from Pegasus.api import *

# --- Available dictionary keys for the OBS fetcher
SOURCES = ['USGS','USGS_RIVERS','NOAA','CONTRAILS','NDBC','NDBC_HISTORIC']
DICT_SOURCES={'USGS_RIVERS':'USGS_RIVERS','USGS':'USGS','NOAA_STATIONS':'NOAAWEB','NDBC_HISTORIC':'NDBC_HISTORIC','NDBC_BUOYS':'NDBC','CONTRAILS_COASTAL':'CONTRAILS','CONTRAILS_RIVERS':'CONTRAILS'}

source_keys={'NOAA':'NOAA_STATIONS',
            'NDBC':'NDBC_BUOYS',
            'NDBC_historic':'NDBC_HISTORIC',
            'CONTRAILS_coastal':'CONTRAILS_COASTAL',
            'CONTRAILS_rivers':'CONTRAILS_RIVERS',
            'USGS_coastal':'USGS',
            'USGS_rivers':'USGS_RIVERS'
            }

# --- Fetch input yaml name and the associated station list filename ( must be local to the yaml file)

def return_inputs_stationnames(yamldir, source_map, source): # : source):
    yaml_file = source_map[source]
    source_yaml=f'{yamldir}/{yaml_file}'
    with open(source_yaml, 'r') as stream:
        config = yaml.safe_load(stream)
        print(f'Opened yaml file {source_yaml}')
    if len(config['SOURCEMAP'].keys()) > 1:
        print(f'Can only use yamls with a single key. You have: {config["SOURCEMAP"].keys()}')
        sys.exit(1)
    basecsv = config['SOURCEMAP'][source_keys[source]]['STATION_FILE'] # Too bad we didn;t simplify the 'NOAA_STATIONS' to just 'STATIONS'
    source_csv = f'{yamldir}/{basecsv}'
    return source_yaml, source_csv

def return_inputs_yaml_stationnames(yamldir, inyaml): # : source):
    source_yaml=inyaml # f'{yamldir}/{yaml_file}'
    source = extract_data_source_name(source_yaml)
    with open(source_yaml, 'r') as stream:
        config = yaml.safe_load(stream)
        print(f'Opened yaml file {source_yaml}')
    if len(config['SOURCEMAP'].keys()) > 1:
        print(f'Can only use yamls with a single key. You have: {config["SOURCEMAP"].keys()}')
        sys.exit(1)
    basecsv = config['SOURCEMAP'][source_keys[source]]['STATION_FILE'] # Too bad we didn;t simplify the 'NOAA_STATIONS' to just 'STATIONS'
    source_csv = f'{yamldir}/{basecsv}'
    return source_yaml, source_csv

# --- Predict the names and numbers of output files. This is tricky because the number of output files depends on how many
#     data products the user specified in the map_source yaml. Also, the filename metadata depends on the stoptime and ndays.
#     In order to properly StgeOut those files we need to know what they all are beforehand.
#     The output filename wll come as PAIR of data and metadata for each product specified
# Eg: noaaweb_stationdata_meta_water_level_2025-08-01T00:00:00.csv  noaaweb_stationdata_water_level_2025-08-01T00:00:00.csv

# TODO refactor this code as it also exists in the run_harverter code but is needed to specify replica and job arguments

def construct_root_filenames(data_source):
    """
    Define (partial) metadata for the current source
    """
    output_fileroot=f'{data_source.lower()}_stationdata'
    output_metafileroot=f'{data_source.lower()}_stationdata_meta'
    return output_fileroot,output_metafileroot

def construct_metadata(data_product, endt):
    """
    Construct remaining file metadata and include the stoptime values needed for DB uploads
    Complexity caused by this provider segregating rivers from coastal data
    """
    if 'river_water_level' in data_product:
         metadata=f'{data_product}_{endt}_RIVERS'
    elif 'river_stream_elevation' in data_product:
         metadata=f'{data_product}_{endt}_RIVERS'
    elif 'coastal_water_level' in data_product:
         metadata=f'{data_product}_{endt}_COASTAL'
    else:
        metadata=f'{data_product}_{endt}'
    return metadata

def works_return_list_outputfilenames(stoptime,source_yaml):
    """
    The order of output filenames is irrelevant as these will all be specified to the Job arguments and order doesn't matter
    """
    dataf = list()
    metaf = list()

    endt=stoptime.replace(' ','T') # For filename metadata
    with open(source_yaml, 'r') as stream:
        source_config = yaml.safe_load(stream)
        data_sources = list(source_config['SOURCEMAP'].keys())
        for data_source_in in data_sources:
            data_source=DICT_SOURCES[data_source_in]
            source_products = list(source_config['SOURCEMAP'][data_source_in]['SOURCES'].values())
            data_source_short = 'CONTRAILS' if 'CONTRAILS' in data_source_in else DICT_SOURCES[data_source_in]
            output_fileroot,output_metafileroot = construct_root_filenames(data_source_short)
            for data_product in source_products:
                metadata = construct_metadata(data_product, endt)
                data_source_short = 'CONTRAILS' if 'CONTRAILS' in data_source_in else data_source
                output_fileroot,output_metafileroot = construct_root_filenames(data_source_short)
                metadata = construct_metadata(data_product, endt)
                print(f'{data_product},{data_source},{data_source_short},{metadata},{output_fileroot},{output_metafileroot}')
                dataf.append(f'{output_fileroot}_{metadata}.csv')
                metaf.append(f'{output_metafileroot}_{metadata}.csv')
    print('Finished with filename')
    return dataf, metaf

def return_list_outputfilenames(stoptime,source_yaml):
    """
    The order of output filenames is irrelevant as these will all be specified to the Job arguments and order doesn't matter
    """
    dataf = list()
    metaf = list()
    bstrings = list()

    endt=stoptime.replace(' ','T') # For filename metadata
    with open(source_yaml, 'r') as stream:
        source_config = yaml.safe_load(stream)
        data_sources = list(source_config['SOURCEMAP'].keys())
        for data_source_in in data_sources:
            data_source=DICT_SOURCES[data_source_in]
            source_products = list(source_config['SOURCEMAP'][data_source_in]['SOURCES'].values())
            #output_fileroot,output_metafileroot = construct_root_filenames(data_source_in)
            data_source_short = 'CONTRAILS' if 'CONTRAILS' in data_source_in else DICT_SOURCES[data_source_in]
            output_fileroot,output_metafileroot = construct_root_filenames(data_source_short)
            for data_product in source_products:
                metadata = construct_metadata(data_product, endt)
                data_source_short = 'CONTRAILS' if 'CONTRAILS' in data_source_in else data_source
                output_fileroot,output_metafileroot = construct_root_filenames(data_source_short)
                metadata = construct_metadata(data_product, endt)
                bstrings.append(f'{output_fileroot}_{data_product}')
                bstrings.append(f'{output_metafileroot}_{data_product}')
                dataf.append(f'{output_fileroot}_{metadata}.csv')
                metaf.append(f'{output_metafileroot}_{metadata}.csv')
    print('Finished with filename')
    return dataf, metaf, bstrings

def return_list_files(indir,ext='yaml'):
    """
    Need to return the list of yamls to be used for processing.
    The plan is for the Pegasus workflow to process all observed *.yamls
    concurrently
    """
    import glob, os
  
    files = list()
    #os.chdir(dir)
    for file in glob.glob(f'{indir}/*.{ext}'):
        files.append(f'{file}')
    return files

def get_key_from_value(input_dict, target_value):
    print(input_dict)
    print(target_value)
    for key, value in input_dict.items():
        print(f' {key}:{value}:{target_value}')
        if value == target_value:
            return key
    return None  

def extract_data_source_name(yamlfile):

    """
    Only return the first value of the list
    We want to get the Key value for the given data_sources for which there should only be one
    """
    try:
        with open(yamlfile, 'r') as stream:
            source_config = yaml.safe_load(stream)
            data_sources = list(source_config['SOURCEMAP'].keys())
            print(source_config)
            print(data_sources)
            if len(data_sources) >1: print('die')
            source = get_key_from_value(source_keys, data_sources[0])
    except Exception as e:
        print(f'Failed to process yaml {yamlfile}: {e}')
    return source

