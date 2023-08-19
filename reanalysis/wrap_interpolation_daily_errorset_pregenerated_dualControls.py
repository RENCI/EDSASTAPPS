#!/usr/bin/env python

import os
import sys
import pandas as pd
import json
from utilities.utilities import utilities
import io_utilities.io_utilities as io_utilities

#############################################################################
# Build a slurm file
#def build_slurm(key, daily_file_errors, main_yamlname, input_directory, output_directory, map_file, gridname, cv_testing):
def build_slurm(key,slrvalue, midvalue, vlfvalue, main_yamlname, input_directory, output_directory, map_file, gridname, cv_testing): 
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
    slurm.append('python -u $dir/interpolation_daily_errorset_pregenerated_dualControls.py --input_directory "'+input_directory+'" --output_directory "'+output_directory+'" --slr_file_errors "'+slrvalue+'" --mid_file_errors "'+midvalue+'" --vlf_file_errors "'+vlfvalue+'"  --main_yamlname "'+main_yamlname+'" --map_file "'+map_file+'" --iometadata "'+key+'" --gridname "'+gridname+'"' )
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
    return ('./tmp/runSlurm'+key+'.sh') # Just use keyvlf as it should be te same for all models

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

    # Pull out the individual model terms from the input dict
    slr_all_dict=files_dict['SLR']
    mid_all_dict=files_dict['MID']
    vlf_all_dict=files_dict['VLF']

    # assemble inputs for a specific run
    # Grab the SLR/MID/VLF and build a date specific new dict for passing to the interpolators.
    # This is looping over DATES and grabs the appropriate filenames to pass to the interpolator

    print(' Performing lOOP xxxxxxxxxxxxxxxxxxxxxxxxxxx')
    print(f' vlfdict {vlf_all_dict}')
    print('TEST ZIP')
    count=0
    for keyslr,keymid,keyvlf in zip(slr_all_dict,mid_all_dict,vlf_all_dict):
        count +=1
        print(f' count {count}')

    for keyslr,keymid,keyvlf in zip(slr_all_dict,mid_all_dict,vlf_all_dict):
        slrvalue=slr_all_dict[keyslr]
        midvalue=mid_all_dict[keymid]
        vlfvalue=vlf_all_dict[keyvlf]
        print(f' keyvlf {keyvlf}')
        print(f' slrval {slrvalue}')
        print(f' mid {keymid}')
        print(f' vlf val {midvalue}')
        print(vlfvalue)
        slurm_filename = build_slurm(keyvlf,slrvalue, midvalue, vlfvalue, main_yamlname, topdir, outputdir, map_file, gridname, cv_testing)
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
                        help='Boolean: Invoke a CV procedure subsequent to fitting interpolation model')
    parser.add_argument('--map_file', action='store',dest='map_file', default=None,
                        help='str: Select appropriate map_file ym; for grid lookup')
    args = parser.parse_args()
    sys.exit(main(args))
