"""
TODO: Instead of finding my old files, put in the actual creation methods here
"""

import pandas as pd

###############################################################################
# run
###############################################################################

def main(out_directory, sv_refactored_path, max_objects=50000000):

    # load up files
    hdf = pd.HDFStore(sv_refactored_path)

    catalog = hdf['/IM3SHAPE/catalog']
    # rename
    catalog.rename(columns={'ra': 'RA', 'dec': 'DEC',
                            'Z_MCMC_0_annz2': 'Z_MCMC_ANNZ2',
                            'Z_MEAN_annz2': 'Z_MEAN_ANNZ2',
                            'Z_MCMC_0_tpz': 'Z_MCMC_TPZ',
                            'Z_MEAN_tpz': 'Z_MEAN_TPZ',
                            'Z_MCMC_0_skynet': 'Z_MCMC_SKYNET',
                            'Z_MEAN_skynet': 'Z_MEAN_SKYNET',
                            'Z_MCMC_0_bpz': 'Z_MCMC_BPZ',
                            'Z_MEAN_bpz': 'Z_MEAN_BPZ',
                            }, inplace=True)

    # drop columns
    catalog = catalog.drop(
        [col for col in catalog.columns if 'ztomo' in col] +
        [col for col in catalog.columns if 'kmeans' in col] +
        [col for col in catalog.columns if '_1_' in col] +
        [col for col in catalog.columns if '_2_' in col] +
        [col for col in catalog.columns if 'MAG_AUTO' in col] +
        ['IM3SHAPE_FLAG', 'NGMIX_FLAG', 'SVA1_FLAG',],
        axis = 1)

    # cut any catalog with photoz_bin == -99
    catalog = catalog[catalog['PHOTOZ_BIN'] != -99].reset_index()

    # save catalog
    catalog.to_hdf(out_directory + 'im3shape_galaxies.h5', 'catalog')

    randoms = hdf['/IM3SHAPE/random']

    randoms.rename(columns={'z': 'Z', 'ra': 'RA', 'dec': 'DEC'}, inplace=True)
    randoms = randoms.drop(
        [col for col in randoms.columns if 'ztomo' in col] +
        [col for col in randoms.columns if 'kmeans' in col],
        axis = 1)

    # assign PHOTOZ_BIN
    randoms['PHOTOZ_BIN'] = -1
    randoms['PHOTOZ_BIN'][(randoms['Z'] > 0.3) & (randoms['Z'] < 0.55)] = 0
    randoms['PHOTOZ_BIN'][(randoms['Z'] > 0.55) & (randoms['Z'] < 0.83)] = 1
    randoms['PHOTOZ_BIN'][(randoms['Z'] > 0.83) & (randoms['Z'] < 1.3)] = 2

    # go from Z to Z_MEAN_BPZ, Z_MCMC_BPZ, etc
    for col1 in ['Z_MEAN_', 'Z_MCMC_']:
        for col2 in ['TPZ', 'ANNZ2', 'SKYNET', 'BPZ']:
            randoms[col1 + col2] = randoms['Z']

    # save randoms
    randoms.to_hdf(out_directory + 'im3shape_randoms.h5', 'catalog')

    ###########################################################################
    # redmapper
    ###########################################################################

    # load redmapper from old catalog
    redmapper = hdf['/REDMAPPER_OLD/catalog']
    # rename
    redmapper.rename(columns={'z': 'Z', 'ra': 'RA', 'dec': 'DEC'}, inplace=True)

    # drop useless stuff
    redmapper = redmapper.drop(
            [col for col in redmapper.columns if 'kmeans' in col]
            + [col for col in redmapper.columns if 'ztomo' in col]
            + [col for col in redmapper.columns if 'zref' in col]
            + ['healpix_region', 'split_region'],
            axis=1)

    # save as hdf
    redmapper.to_hdf(out_directory + 'redmapper_galaxies.h5', 'catalog')

    # randoms
    redmapper_randoms = hdf['/REDMAPPER_OLD/random']
    # rename
    redmapper_randoms.rename(columns={'z': 'Z', 'ra': 'RA', 'dec': 'DEC'}, inplace=True)
    redmapper_randoms = redmapper_randoms.drop(
            [col for col in redmapper_randoms.columns if 'kmeans' in col]
            + [col for col in redmapper_randoms.columns if 'ztomo' in col]
            + [col for col in redmapper_randoms.columns if 'zref' in col]
            + ['healpix_region', 'split_region'],
            axis=1)

    redmapper_randoms.to_hdf(out_directory + 'redmapper_randoms.h5', 'catalog')

if __name__ == '__main__':
    out_directory = '.'
    refactor_sv_path = 'data.h5'
    main(out_directory, refactor_sv_path)
