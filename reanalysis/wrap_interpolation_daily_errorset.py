#!/usr/bin/env python

import os
import sys
import pandas as pd
import json
from utilities.utilities import utilities
import io_utilities.io_utilities as io_utilities

#############################################################################
# Build a slurm file

# slurm.append('python -u $dir/krigListOfErrorSets.py --daily --inrange "'+RANGE+'" --insill "'+SILL+'" --outroot "'+ROOTDIR+'" --yamlname "'+YAMLNAME+'" --errorfile "'+ERRFILE+'" --clampfile "'+CLAMPFILE+'" --controlfile "'+CONTROLFILE+'" --gridjsonfile "'+ADCJSON+'"' )

def build_slurm(key, daily_file_errors, main_yamlname, input_directory, output_directory, map_file, gridname, cv_testing): 
    slurm = list()
    slurm.append('#!/bin/sh')
    slurm.append('#SBATCH -t 24:00:00')
    slurm.append('#SBATCH -p batch')
    slurm.append('#SBATCH -N 1')
    slurm.append('#SBATCH -n 1')
    slurm.append('#SBATCH -J Interpolate')
    slurm.append('#SBATCH -p lowpri')
    slurm.append('#SBATCH --mem-per-cpu 64000')
    slurm.append('echo "Begin the Interpolation phase" ')
    slurm.append('export PYTHONPATH=/projects/sequence_analysis/vol1/prediction_work/AST:/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS')
    slurm.append('dir="/projects/sequence_analysis/vol1/prediction_work/EDSASTAPPS/reanalysis"')
    slurm.append('python -u $dir/interpolation_daily_errorset.py --cv_testing --input_directory "'+input_directory+'" --output_directory "'+output_directory+'" --daily_file_errors "'+daily_file_errors+'" --main_yamlname "'+main_yamlname+'" --map_file "'+map_file+'" --iometadata "'+key+'" --gridname "'+gridname+'"' )
    try:
        os.makedirs('./tmp')
    except OSError as exc: # Python >2.5
        #if exc.errno == errno.EEXIST and os.path.isdir('./tmp'):
        pass
        #else: raise
    with open('./tmp/runSlurm'+key+'.sh', 'w') as file:
        for row in slurm:
            file.write(row+'\n')
    file.close()
    return ('./tmp/runSlurm'+key+'.sh')

##
## Read the relevant input json file and loop over the contents
##

def main(args):
    #assumes the existance of a proper main.yml to get IO information
    if args.main_yamlname is None:
        main_yamlname=os.path.join(os.path.dirname(__file__), './config', 'main.yml')
    else:
        main_yamlname=args.main_yamlname
    config = utilities.init_logging(subdir=None, config_file=main_yamlname)
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
    map_file=args.map_file
    
    filepath = args.json_daily_file_errors
    files_dict = io_utilities.read_json_file(filepath) 
    utilities.log.info('List of error files found {}'.format(files_dict))

    # assemble inputs for a specific run

    for key, value in files_dict.items():
        print(key)
        print(value)
        print('Start {}'.format(key))
        slurm_filename = build_slurm(key, value, main_yamlname, topdir, outputdir, map_file, gridname, cv_testing)
        cmd = 'sbatch ./'+slurm_filename
        print(cmd)
        os.system(cmd)

    print('Completed ensemble')

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--main_yamlname', action='store',dest='main_yamlname', default=None,
                        help='str: Select appropriate main_yamlname')
    parser.add_argument('--json_daily_file_errors', action='store', dest='json_daily_file_errors', default=None,
                        help='JSON that contains the complete set of compute_error files')
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
    args = parser.parse_args()
    sys.exit(main(args))
