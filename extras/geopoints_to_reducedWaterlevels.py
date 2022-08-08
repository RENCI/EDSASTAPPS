#!/usr/bin/env python
# coding: utf-8

# Sample URLs at:
#url='/projects/sequence_analysis/vol1/prediction_work/NOAA.Reanalysis/DATATEST/BIG_TEST/fort.63_transposed_and_rechunked_1024.nc'
#url='/projects/reanalysis/ADCIRC/ERA5/hsofs/2019/fort.63_transposed_and_rechunked_1024.nc'
#url='/projects/sequence_analysis/vol1/prediction_work/NOAA.Reanalysis/DATATEST/BIG_TEST/swan_HS.63_2019_transposed_and_rechunked_1024.nc'


import sys
import pandas as pd
import numpy as np
import time as tm
import adcirc_utilities as utilities

def main(args):

    url = args.url
    variable_name=args.variable_name
    geopointsfile=args.geopointsfile

    if url is None:
        url='/projects/reanalysis/ADCIRC/ERA5/hsofs/2019/fort.63_transposed_and_rechunked_1024.nc'
        variable_name='zeta'
        print('Selecting a test data set')

    if args.geopointsfile is None:
        geopointsfile='./testdata/hsofs_example_geopoints.csv'

    nearest_neighbors=args.kmax

    dump_header=False

    print(f'input URL selected {url}')
    print(f'Indicated var name is {variable_name}')
    print(f'Geopoints coming from the file {geopointsfile}')
    print(f'Selected nearest neighbors values is {nearest_neighbors}')
    print(f'Status to dump headers is {dump_header}')

##
## Read the geopoints data and extract the lon/lat data
##

    #df_geopoints = pd.read_csv(geopointsfile, index_col=0, header=0, skiprows=[1])
    df_geopoints = pd.read_csv(geopointsfile, index_col=0, header=0)

    print(df_geopoints) # Keep this for later comparisons.
    geopoints = df_geopoints[['lon','lat']].to_numpy()
    print(f'Number of input geopoints is {len(geopoints)}')

##
## Initialize input the ADCIRC and popuate with grid properties
##

    ds=utilities.f63_to_xr(url)
    agdict=utilities.get_adcirc_grid_from_ds(ds)
    agdict=utilities.attach_element_areas(agdict)

##
## Perform the Tree building, Query and Weighted averaging of the timeseries (reduction)
##

    agdict=utilities.ComputeTree(agdict)
    agresults=utilities.ComputeQuery(geopoints, agdict, kmax=nearest_neighbors)
    agresults=utilities.ComputeBasisRepresentation(geopoints, agdict, agresults)
    agresults=utilities.ConstructReducedWaterLevelData_from_ds(ds, agdict, agresults, variable_name=variable_name)

##
## For future ref what elements where not assigned to any element
##

    print(f'List of {len(agresults["outside_elements"])} stations not assigned to any grid element follows for kmax {nearest_neighbors}')
    print(df_geopoints.iloc[agresults['outside_elements']])

##
## Assemble data for storage
##

    df_product_data=agresults['final_reduced_data']
    df_product_metadata=agresults['final_meta_data']
    df_product_data.to_csv('data.csv',header=True)
    df_product_metadata.to_csv('meta.csv',header=True)

    print(df_product_data)
    print('Finished')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--url', action='store', dest='url', default=None, type=str,
                        help='url to fetch ADCIRC netCDF data')
    parser.add_argument('--geopointsfile', action='store', dest='geopointsfile', default=None, type=str,
                        help='filename to find geopoints lons/lats file data')
    parser.add_argument('--variable_name', action='store', dest='variable_name', default=None, type=str,
                        help='Variable name of interest from the supplied url')
    parser.add_argument('--kmax', action='store', dest='kmax', default=10, type=int,
                        help='nearest_neighbors values when performing the Query')
    args = parser.parse_args()
    sys.exit(main(args))
