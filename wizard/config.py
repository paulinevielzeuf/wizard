# TODO: License

"""
.. module:: config
"""

from __future__ import print_function
from os import path, makedirs
import numpy as np

from .dataset import dataset
from .paircounts import paircounts
from .combine import combine
from .dndz import dndz
from .plot import plot
from .compare_noz_fom import compare_noz_fom
from .FoM import FoM
from .plot_FoM import plot_FoM


def read_config(file_name):
    """Read a configuration dictionary from a file

    :param file_name:   yaml file name which we read
    """
    import yaml
    with open(file_name) as f_in:
        config = yaml.load(f_in.read())
    return config

def save_config(config, file_name):
    """Take a configuration dictionary and save to file

    :param config:  Dictionary of configuration
    :param file_name: file we are saving to
    """
    import yaml
    with open(file_name, 'w') as f:
        f.write(yaml.dump(config, default_flow_style=False))

def make_directories(config):
    # go through the known dataset paths and make directories
    paths = {'dataset': ['dataset_path'],
             'paircounts': ['paircounts_path'],
             'combine': ['combine_path'],
             'dndz': ['dndz_path'],
             'plot': ['plot_path'],
             'compare': ['compare_path'],
             'compare_noz_fom': ['compare_noz_fom_path'],
              'FoM' : ['FoM_path'],
              'plot_FoM' : ['plot_path'],}
    for i in paths:
        if i in config:
            for j in paths[i]:
                if j in config[i]:
                    check_make(path.dirname(config[i][j]))

def check_make(path_check):
    """
    Convenience routine to avoid that annoying 'can't make directory; already
    present!' error.
    """
    if not path.exists(path_check):
        makedirs(path_check)

def wizard(config):
    """Run the wizard code, given a configuration dictionary

    Parameters
    ----------
    config : the dictionary of config files

    """
    from time import time
    if config['verbose']:
        time0 = time()
    else:
        time0 = 0

    for entry in config['run']:
        if entry == 'make_directories':
            if time0:
                print('make directories', time() - time0)
            make_directories(config)
        elif entry == 'dataset':
            if time0:
                print('dataset', time() - time0)
                if config['verbose'] > 1:
                    time0i = time()
                else:
                    time0i = 0
            reference_catalogs, unknown_catalogs = dataset(time0=time0i, **config['dataset'])
        elif entry == 'paircounts':
            """
            We need from datasets:
            reference_catalogs [output form dataset function]
            unknown_catalogs [output from dataset function]
            reference_bins [output from config]
            unknown_bins [output from config]
            """
            if time0:
                print('paircounts', time() - time0)
                if config['verbose'] > 1:
                    time0i = time()
                else:
                    time0i = 0
            pairs = paircounts(reference_catalogs, config['dataset']['reference_bins'],
                               unknown_catalogs, config['dataset']['unknown_bins'],
                               config['dataset']['healpix_nside'], time0=time0i, **config['paircounts'])
        elif entry == 'combine':
            """
            pairs and healpix_nside come from previous run of paircounts
            """
            if time0:
                print('combine', time() - time0)
                if config['verbose'] > 1:
                    time0i = time()
                else:
                    time0i = 0
            if 'paircounts' not in config['run']:
                # need to come up with pairs as a default
                pairs = {'DuDr': None, 'RuRr': None, 'DuRr': None, 'RuDr': None}
            combined_pairs = combine(reference_catalogs, config['dataset']['reference_bins'], unknown_catalogs, config['dataset']['unknown_bins'], pairs, config['dataset']['healpix_nside'], time0=time0i, **config['combine'])
        elif entry == 'dndz':
            if time0:
                print('dndz', time() - time0)
                if config['verbose'] > 1:
                    time0i = time()
                else:
                    time0i = 0
            if 'combine' not in config['run']:
                # put in false combined_pairs since it probably doesn't matter
                combined_pairs = None
            min_sep = config['paircounts']['min_sep']
            max_sep = config['paircounts']['max_sep']
            nbins = config['paircounts']['nbins']
            redshift_bins = np.arange(*config['dataset']['reference_bins'][config['dndz']['redshift_column']][1])
            redshifts = np.append(-10, np.append(0.5 * (redshift_bins[1:] + redshift_bins[:-1]), 10))
            phi, cov, stats, tags = dndz(combined_pairs, redshifts, min_sep=min_sep, max_sep=max_sep, nbins=nbins, time0=time0i, **config['dndz'])
        elif entry == 'plot':
            # plot the results
            if time0:
                print('plot', time() - time0)
                if config['verbose'] > 1:
                    time0i = time()
                else:
                    time0i = 0
            z_min = config['dndz']['z_min']
            z_max = config['dndz']['z_max']


            # redshift_column = config['dndz']['redshift_column']
            plot(phi, cov, stats, tags, unknown_catalogs[0], redshift_bins,
                 z_min=z_min, z_max=z_max, #redshift_column=redshift_column,
                 time0=time0,
                 **config['plot'])

        elif entry == 'compare':

            # fit to another catalog

            # repeat plotting
            pass
        elif entry == 'compare_noz_fom':
                if time0:
                    print('compare_noz_fom', time() - time0)
                    if config['verbose'] > 1:
                            time0i = time()
                    else:
                            time0i = 0
                z_min = config['dndz']['z_min']
                z_max = config['dndz']['z_max']
                path_rec=config['dndz']['dndz_path']
                biaised_nz=compare_noz_fom(phi,tags,unknown_catalogs[0], redshift_bins,
                z_min=z_min, z_max=z_max, path_rec=path_rec,#redshift_column=redshift_column,
                time0=time0,
                **config['compare_noz_fom'])
        elif entry == 'FoM':
                if time0:
                    print('compare_noz_fom', time() - time0)
                    if config['verbose'] > 1:
                            time0i = time()
                    else:
                            time0i = 0
                cosmosissource = 'source my-source'
                Nz_path=config['compare_noz_fom']['compare_noz_fom_path']
                fom=FoM(Nz_path=Nz_path,time0=time0,**config['FoM'])
        elif entry == 'plot_FoM':
                if time0:
                    print('compare_noz_fom', time() - time0)
                    if config['verbose'] > 1:
                            time0i = time()
                    else:
                            time0i = 0
                fom_path=config['FoM']['FoM_path']
                comp=config['compare_noz_fom']['redshift_catalog_columns']
                bins=config['dataset']['unknown_bins']
                plot_fom=plot_FoM(fom_path=fom_path,comp=comp,bins=bins,time0=time0,**config['plot_FoM'])

    if time0:
        print('done', time() - time0)
