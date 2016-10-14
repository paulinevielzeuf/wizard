"""
Point of this script is to submit a ton of batch jobs for the healpixels.
After all, why use only one computer when you have several thousand at your
fingertips?

Note that you should first run with --test command, which will run one healpix
iteration after making the dataset. This lets you determine the amount of
memory and time you should request to your batch farm!
"""
from __future__ import division, print_function

from subprocess import call, check_output
import pandas as pd
from math import ceil

from wizard import save_config, read_config, check_make, wizard

config = read_config('config.yaml')

###############################################################################
# config parameters -- ADJUST
###############################################################################

config_name = 'refREDMAPPERSV_unkIM3SHAPE'
sv_refactored_path = '/nfs/slac/g/ki/ki18/des/cpd/cluster-z/refactor_sv/data.h5'
make_dataset_kwargs = {'sv_refactored_path': sv_refactored_path}

# batch farm parameters
memory = 5000 # Mb
req_time_per_pixel = 100 # minutes PER healpixel
only_run = 200  # maximum number of concurrent jobs

# batch farm command generator
job_directory = '/nfs/slac/g/ki/ki19/des/cpd/wizard/{0}/'.format(config_name)
def batch_commander(label, configfile, maxnum, only_run=200, check=True):
    # batch farm memory requirements
    req_time = req_time_per_pixel * config['dataset']['healpix_filter_size']

    logfile = '{0}/logs/{1}.log'.format(job_directory, 'paircounts')

    jobname = '"wizard_{0}[1-{1}]%{2}"'.format(config_name, maxnum, only_run)

    command = ['bsub',
               '-J', jobname,
               '-o', logfile,
               '-W', str(req_time),
               '-M', str(memory),
               '-R', '"span[hosts=1] rusage[mem={0}]"'.format(memory),
               'wizard', configfile]

    # check if we want to run this job
    if check:
        jobcheck = check_output(['bjobs', '-wJ', jobname])
        if jobname in jobcheck:
            command = ['echo', 'skipping {0} because it is already running'.format(label)]

    return command

###############################################################################
# Update config to include job_directory
###############################################################################

config['verbose'] = 2

config['dataset']['reference_data_path'][0] = job_directory + config['dataset']['reference_data_path'][0]
config['dataset']['reference_random_path'][0] = job_directory + config['dataset']['reference_random_path'][0]
config['dataset']['unknown_data_path'][0] = job_directory + config['dataset']['unknown_data_path'][0]
config['dataset']['unknown_random_path'][0] = job_directory + config['dataset']['unknown_random_path'][0]
config['dataset']['dataset_path'] = job_directory + config['dataset']['dataset_path']

config['paircounts']['paircounts_path'] = job_directory + config['paircounts']['paircounts_path']
config['combine']['combine_path'] = job_directory + config['combine']['combine_path']
config['dndz']['dndz_path'] = job_directory + config['dndz']['dndz_path']
config['plot']['plot_path'] = job_directory + config['plot']['plot_path']
config['compare']['compare_path'] = job_directory + config['compare']['compare_path']

###############################################################################
# check if we are in testing mode
###############################################################################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--test', action='store_true', dest='test', help='Do a single healpix')
parser.add_argument('--call', action='store_true', dest='call', help='Submit healpixels to batch farm')
parser.add_argument('--dataset', action='store_true', dest='dataset', help='Create the dataset')
args = parser.parse_args()

###############################################################################
# make some directories
###############################################################################

# need to make a dataset directory where the new datasets will live
check_make('{0}/dataset'.format(job_directory))
# need to make a config directory where the oodles of new configs live
check_make('{0}/configs'.format(job_directory))
# need to make a logs directory to check out how my jobs performed
check_make('{0}/logs'.format(job_directory))

###############################################################################
# make fits dataset
###############################################################################

if args.dataset:
    from make_dataset import main
    main(job_directory + '/dataset/', **make_dataset_kwargs)

###############################################################################
# make wizard dataset (needed to determine number of jobs to submit)
###############################################################################

# edit config to only run make_directories
make_directories_config = config.copy()
make_directories_config['run'] = ['make_directories']
# execute config file with wizard
wizard(make_directories_config)

if args.dataset:
    # edit config to only run dataset
    dataset_config = config.copy()
    dataset_config['run'] = ['make_directories', 'dataset']
    # execute config file with wizard
    wizard(dataset_config)

###############################################################################
# create and execute config files
###############################################################################

# first load up the dataset to get a sense of the number of healpixels involved
hdf = pd.HDFStore(config['dataset']['dataset_path'])
# now it _should_ be that the randoms and data have the same healpixels in this
# coarse limit... but let's go with randoms for now
data = hdf['/reference/random']
hpix = data['HPIX'].unique()

# now iterate
# check if we are in testing mode
if args.test:
    ith_iter = [0]
else:
    ith_iter = range(int(ceil(len(hpix) / config['dataset']['healpix_filter_size'])))

for ith in ith_iter:
    # adjust config
    config_run = config.copy()
    config_run['run'] = ['dataset', 'paircounts']
    config_run['dataset']['load_dataset'] = True
    config_run['dataset']['healpix_filter_start'] = ith

    # write config
    configfile = '{0}/configs/{1}.yaml'.format(job_directory, ith + 1)
    save_config(config_run, configfile)

# get config command
command = batch_commander('\$LSB_JOBINDEX', '{0}/configs/\$LSB_JOBINDEX.yaml'.format(job_directory), maxnum=max(ith_iter), check=False)

# execute config
print(' '.join(command))
if args.call:
    # TODO: so, uh, it will not work if you call from in here but if you copy and paste the command it works?!
    call(command)

###############################################################################
# write out config file for combining
###############################################################################

config_combine = config.copy()
config_combine['run'] = ['dataset', 'paircounts', 'combine']
config_combine['dataset']['load_dataset'] = True
config_combine['dataset']['healpix_filter_start'] = -1
config_combine['paircounts']['load_paircounts'] = True
configfile = '{0}/configs/combine.yaml'.format(job_directory)
save_config(config_combine, configfile)

###############################################################################
# write out config file for just dndz
###############################################################################

config_dndz = config.copy()
config_dndz['run'] = ['dataset', 'combine', 'dndz']
# speed up dataset by getting rid of paths for randoms
config_dndz['dataset']['reference_random_path'] = ''
config_dndz['dataset']['unknown_random_path'] = ''
config_dndz['dataset']['load_dataset'] = True
config_dndz['dataset']['healpix_filter_start'] = -1
config_dndz['paircounts']['load_paircounts'] = True
config_dndz['combine']['load_combine'] = True
configfile = '{0}/configs/dndz.yaml'.format(job_directory)
save_config(config_dndz, configfile)

###############################################################################
# write out config file for just analysis
###############################################################################

config_analyze = config.copy()
config_analyze['run'] = ['dataset', 'dndz', 'plot', 'compare']
# speed up dataset by getting rid of paths for randoms
config_analyze['dataset']['reference_random_path'] = ''
config_analyze['dataset']['unknown_random_path'] = ''
config_analyze['dataset']['load_dataset'] = True
config_analyze['dataset']['healpix_filter_start'] = -1
config_analyze['paircounts']['load_paircounts'] = True
config_analyze['combine']['load_combine'] = True
config_analyze['dndz']['load_dndz'] = True
configfile = '{0}/configs/analyze.yaml'.format(job_directory)
save_config(config_analyze, configfile)
