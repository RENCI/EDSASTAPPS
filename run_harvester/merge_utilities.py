#!/usr/bin/env python

##############
## Merge one or more data files of a specific SOURCE into a single dataframe of dims (ntimes,nstations)
## On entry, we must have already determined a proper metadata string (bstring) to search out the proper list of files
## for processing
##############

import os, sys
import pandas as pd
from pathlib import Path
from utilities.utilities import utilities as utilities

class merge_utilities(object):
    obs_filelist = None
    outfilename = None

    def __init__(self, obs_filelist, outfilename):
        """
        """
        self.obs_filelist = obs_filelist
        self.outfilename = outfilename

    def aggregate_source_dataproduct(self, filelist):
        """
        We expect only a single data_product and the output filename has already been fully specified. No more bstrings.
        """
        df_list = list()
        #for bfile in [name for name in filelist if bstring in name]:
        for bfile in filelist:
            df = pd.read_csv(bfile, header=0, index_col=0)
            if not 'meta' in bfile:
                df.index = pd.to_datetime(df.index)
            df_list.append(df)
        df_full = pd.concat(df_list,axis=0)
        df_full = df_full.drop_duplicates(keep='first')
        return df_full

def main(args):
    """
    Aggregate a set of data or metadata files. We expect these data sets to
    be for the same data-product and source
    """
    main_config = utilities.init_logging(subdir=None,config_file='./main.yml')
    outfilename = args.outfilename
    obs_filelist = args.obs_filelist
    utilities.log.info(f'Begin Merger: Input filelist {outfilename}, {args.obs_filelist}')
    utilities.log.info("Product Level Working in {}.".format(os.getcwd()))
    obs_merge = merge_utilities(obs_filelist, outfilename)
    df_full = obs_merge.aggregate_source_dataproduct(obs_filelist)
    df_full.to_csv(outfilename)
    utilities.log.info(f'Merge: Finished writing {outfilename}')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--obs_filelist', nargs="*", required=True, action='store', dest='obs_filelist',
                        help='List of input data or metadata filenames for aggregation'),
    parser.add_argument('--outfilename', required=True, action='store', dest='outfilename', type=str,
                        help='Filename into which merged data will be locally stored')

    args = parser.parse_args()
    sys.exit(main(args))
